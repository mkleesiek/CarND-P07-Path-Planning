#include <iostream>
#include "planner.h"
#include "helpers.h"
#include "spline.h"

using namespace std;

void Planner::reset()
{
  m_current_state = CarState{};
  m_other_cars.clear();
  m_planned_waypoints.clear();
  m_planned_path.clear();
}

void Planner::setMapData(const vector<double>& map_x, const vector<double>& map_y, const vector<double>& map_s,
  double max_s)
{
  // distance of interpolated map waypoints
  constexpr double max_waypoint_distance = 2.0;

  m_map = Map{};

  // drop maps with less than 5 points
  const size_t n_points = min({map_x.size(), map_y.size(), map_s.size()});
  if (n_points < 5) {
    return;
  }

  // if the provided map has a good precision, keep it as is.
  if (max_s / n_points < max_waypoint_distance) {
    m_map = Map{map_x, map_y, map_s, max_s};
    return;
  }

  // otherwise, use C-splines to interpolate a more fine-grained map
  m_map.max_s = max_s;

  vector<double> x_vals, y_vals, s_vals;
  for (int k = -1; k <= 1; k++) {
    for (size_t i = 0; i < n_points; i++) {
      x_vals.push_back(map_x[i]);
      y_vals.push_back(map_y[i]);
      s_vals.push_back(map_s[i] + max_s * k);
    }
  }
  tk::spline sx(s_vals, x_vals);
  tk::spline sy(s_vals, y_vals);

  size_t n_interp_points = max_s / max_waypoint_distance;
  for (size_t i = 0; i < n_interp_points; i++) {
    const double s = max_s * i / n_interp_points;
    m_map.waypoints_s.push_back(s);
    m_map.waypoints_x.push_back(sx(s));
    m_map.waypoints_y.push_back(sy(s));
  }
}

void Planner::setSpeedLimitMPH(double mph)
{
  m_speedLimit = mph_to_ms(mph);
  cout << "Setting speed limit to " << mph << " miles per hour / " << m_speedLimit << " meters per second." << endl;
}

void Planner::updateCarLocalization(double x, double y, double s, double d, double yaw, double speedMPH)
{
  m_current_state.s = s;
  m_current_state.d = d;

  // current time coordinate is unknown, unless the planner has a recent path planned
  // find the point from the previous planned path, which is closest to the current car position
  double min_distance = numeric_limits<double>::infinity();
  auto min_it = m_planned_path.cend();
  for(auto it = m_planned_path.cbegin(); it != m_planned_path.cend(); it++) {
    const double distance = distance_squared(it->x, it->y, x, y);
    if (distance < min_distance) {
      min_distance = distance;
      min_it = it;
    }
  }
  m_current_state.t = (min_it == m_planned_path.cend()) ? 0.0 : min_it->t;

  // estimate s_dot and d_dot from total speed and yaw
  const double speed = mph_to_ms(speedMPH);

  if (fabs(speed) > 0.5) {
    const double x_dot = speed * cos(deg2rad(yaw));
    const double y_dot = speed * sin(deg2rad(yaw));

    const auto sd_2 = getFrenet(x + x_dot * 1.0, y + y_dot * 1.0, deg2rad(yaw));
    m_current_state.s_dot = sd_2.first - s / 1.0;
    m_current_state.d_dot = sd_2.second - d / 1.0;
  }
}

void Planner::updateSensorFusion(const vector<vector<double>>& sensor_data)
{
  m_other_cars.clear();

  // store coordinates and velocities of all other cars
  for (const auto& dataSet : sensor_data)
  {
    const double s = dataSet[5];
    const double d = dataSet[6];

    // approximate frenet velocities from heading and cartesian velocities
    const double heading = atan2(dataSet[4], dataSet[3]);
    const auto sd_2 = getFrenet(dataSet[1] + dataSet[3] * 1.0, dataSet[2] + dataSet[4] * 1.0, heading);
    const double s_dot = sd_2.first - s / 1.0;
    const double d_dot = sd_2.second - d / 1.0;

    m_other_cars.push_back(CarState{s, d, s_dot, d_dot, m_current_state.t});
  }
}

void Planner::planTrajectory(vector<double> &next_x, vector<double> &next_y)
{
  // 1st waypoints for a clean start
  if (m_planned_waypoints.empty() || essentially_equal(m_current_state.t, 0.0, 1E-4)) {
    m_planned_waypoints.clear();
    m_current_state.t = 0.0;
    m_planned_waypoints.push_back(m_current_state);
  }

  // remove previous waypoints which have past
  while (m_planned_waypoints.size() > 3 && m_planned_waypoints[3].t < m_current_state.t - m_buffer) {
    m_planned_waypoints.pop_front();
  }

  // remove previous waypoints lying too far in the future
  while (m_planned_waypoints.size() > 2 && m_planned_waypoints[m_planned_waypoints.size()-2].t > m_current_state.t + m_buffer) {
    m_planned_waypoints.pop_back();
  }

  // keep suggesting and adding new waypoints until the lookahead buffer time is filled
  while (m_planned_waypoints.size() < 5 || m_planned_waypoints.back().t < m_current_state.t + m_buffer + m_look_ahead) {
    CarState next_state = suggestNextWaypoint(m_planned_waypoints.back());
    m_planned_waypoints.push_back(next_state);
  }

  // generate an interpolated jerk-free trajectory from the suggested waypoints
  auto planned_jerkfree = generateJerkMinimizingTrajectory(m_planned_waypoints);

  // collect x/y points for constructing a C-spline
  vector<double> t_vals, s_vals, x_vals, y_vals;

  for (auto& state : planned_jerkfree) {
    const auto xy = getXY(state.s, state.d);
    t_vals.push_back(state.t);
    s_vals.push_back(state.s);
    x_vals.push_back(xy.first);
    y_vals.push_back(xy.second);
  }

  // create C-splines for x and y coordinates
  // we are using C^1 Hermite splines, since updates to the path 2 segments away have no global impact
  tk::spline sx(t_vals, x_vals, tk::spline::cspline_hermite);
  tk::spline sy(t_vals, y_vals, tk::spline::cspline_hermite);

  // use the splines to construct the actual x/y path in 0.02s increments,
  // return those values and store a copy in 'm_planned_path'.
  m_planned_path.clear();

  const double planned_duration = t_vals.back() - t_vals.front();
  const size_t n_interpolated = planned_duration / m_update_interval;

  for (long i = 0; i <= n_interpolated; i++) {

    const double t = t_vals.front() + planned_duration * i / n_interpolated;

    // skip values lying before the current car state
    if (t < m_current_state.t + m_update_interval * 0.5) {
      continue;
    }
    // skip values lying too far in the future
    else if (t > m_current_state.t + m_buffer + m_look_ahead) {
      break;
    }

    PathPoint p{sx(t), sy(t), t};
    next_x.push_back(p.x);
    next_y.push_back(p.y);
    m_planned_path.push_back(std::move(p));
  }
}

vector<Planner::CarState> Planner::generateJerkMinimizingTrajectory(const deque<CarState> &waypoints) const
{
  vector<CarState> jerkmin_traj;
  jerkmin_traj.push_back(waypoints.front());

  for (size_t i = 1; i < waypoints.size(); i++) {

    const auto& last_state = waypoints[i-1];
    auto next_state = waypoints[i];

    // determine start and end conditions of a trajectory segment
    // we assume zero acceleration at the start/end of a segment
    const double traj_duration = next_state.t - last_state.t;
    const auto coeff_s = jmt({last_state.s, last_state.s_dot, 0.0}, {next_state.s, next_state.s_dot, 0.0}, traj_duration);
    const auto coeff_d = jmt({last_state.d, last_state.d_dot, 0.0}, {next_state.d, next_state.d_dot, 0.0}, traj_duration);

    // set the number of intermediate waypoints to generate
    const size_t jmt_iterations = max<size_t>(2, traj_duration / (m_update_interval * 10.0));

    for (size_t k = 0; k < jmt_iterations; k++) {
      const double dt = traj_duration * (k + 1) / jmt_iterations;
      CarState inter_state{
        polynom(dt, coeff_s),
        polynom(dt, coeff_d),
        polynom_gradient(dt, coeff_s),
        polynom_gradient(dt, coeff_d),
        last_state.t + dt
      };
      jerkmin_traj.push_back(inter_state);
    }
  }

  return jerkmin_traj;
}

Planner::CarState Planner::suggestNextWaypoint(const CarState& previous) const
{
  const int current_lane = getLane(previous.d);

  // if we are starting, accelerate straight for 5 seconds
  if (previous.t < 1.0) {
    CarState next_start{};
    next_start.t = 5.0;
    next_start.s_dot = min(5.0 * next_start.t, m_speedLimit);
    next_start.s = previous.s + next_start.s_dot / 2.0 * next_start.t;
    next_start.d = getD(current_lane);
    return next_start;
  }

  // 0) prepare trajectories for going straight and changing to the left/right lane

  CarState next_straight{};
  next_straight.t = previous.t + m_dt_straight;
  next_straight.s_dot = min(previous.s_dot + m_max_accel * m_dt_straight, m_speedLimit);
  next_straight.s = previous.s + (next_straight.s_dot + previous.s_dot) / 2.0 * m_dt_straight;
  next_straight.d = getD(current_lane);
  auto collision_straight = checkCollision(previous, next_straight);

  CarState next_left{};
  next_left.t += previous.t + m_dt_lanechange;
  // reduce target speed and maximum long. acceleration during lane change
  next_left.s_dot = min(previous.s_dot + (m_max_accel * 0.5) * m_dt_lanechange, m_speedLimit * 0.95);
  next_left.s = previous.s + (next_left.s_dot + previous.s_dot) / 2.0 * m_dt_lanechange;
  next_left.d = getD(current_lane - 1);
  auto collision_left = checkCollision(previous, next_left);

  CarState next_right{next_left};
  next_right.d = getD(current_lane + 1);
  auto collision_right = checkCollision(previous, next_right);

  // Evaluate multiple trajectories:
  // 1) if we are on an outer lane, consider changing

  if (current_lane == m_lanes-1 && !collision_left.first) {
    return next_left;
  }
  else if (current_lane == 0 && !collision_right.first) {
    return next_right;
  }

  // 2) keep current lane otherwise (so collision, no slowing down)

  if (!collision_straight.first) {
    return next_straight;
  }

  // 3) if there is a vehicle ahead on the same lane we'd have to slow down in order not to collide,
  //    so consider changing lanes

  // check left-hand lane
  if (!collision_left.first) {
    return next_left;
  }
  // check right-hand lane
  else if (!collision_right.first) {
    return next_right;
  }

  // 4) we can't change lanes right now and have to reduce speed to avoid collision

  // adjust to the speed of the car in front
  next_straight.s_dot = min(next_straight.s_dot, collision_straight.second * 0.9);
  // but don't exceed the (negative) acceleration limit
  next_straight.s_dot = max(previous.s_dot - m_max_accel * m_dt_straight, next_straight.s_dot);
  // estimate distance based on target velocity
  next_straight.s = previous.s + (next_straight.s_dot + previous.s_dot) / 2.0 * m_dt_straight;

  return next_straight;
}

pair<bool, double> Planner::checkCollision(const CarState& me_current, const CarState& me_next) const
{
  const int my_lane_now = getLane(me_current.d);
  const int my_lane_next = getLane(me_next.d);

  // return true for collision if we are attempting to go outside a valid lane
  if (my_lane_next < 0 || my_lane_next >= m_lanes) {
    return {true, 0.0};
  }

  double min_distance = m_distance_to_others;
  CarState candidate{};

  for (const auto& other : m_other_cars) {

    const int others_lane = getLane(other.d);

    // skip vehicles which will not be on my lane
    if (my_lane_next != others_lane) {
      continue;
    }

    // skip other vehicles which are right behind me on my lane!
    if (others_lane == my_lane_now && me_current.t >= other.t && distance_s(me_current.s, other.s) > 0.0) {
      continue;
    }

    // extrapolate the other car position to my next one
    double s_other = other.s + other.s_dot * (me_next.t - other.t);
    double distance_next = distance_s(s_other, me_next.s);

    // check for collision on my next state (being more tolerant on cars behind me)
    if ((distance_next >= 0.0 && distance_next < min_distance)
      || (distance_next < 0.0 && distance_next * -2.0 < min_distance)) {
      candidate = other;
      min_distance = fabs(distance_next);
    }
  }

  return {(min_distance < m_distance_to_others), candidate.s_dot};
}

pair<double, double> Planner::getFrenet(double x, double y, double theta) const
{
  const auto sd = ::getFrenet(x, y, theta, m_map.waypoints_x, m_map.waypoints_y);
  return make_pair(sd[0], sd[1]);
}

pair<double, double> Planner::getXY(double s, double d) const
{
  s = normalize_s(s);
  const auto xy = ::getXY(s, d, m_map.waypoints_s, m_map.waypoints_x, m_map.waypoints_y);
  return make_pair(xy[0], xy[1]);
}

double Planner::normalize_s(double s) const
{
  return fmod(s, m_map.max_s);
}

double Planner::distance_s(double s_front, double s_back) const
{
  s_front = normalize_s(s_front);
  s_back = normalize_s(s_back);

  // increase the front position by one loop length, if necesary
  if (s_front < m_map.max_s * 0.25 && s_back > m_map.max_s * 0.75) {
    s_front += m_map.max_s;
  }

  return s_front - s_back;
}

int Planner::getLane(double d) const
{
  return (d >= 0.0) ? d / m_lane_width : d / m_lane_width - 1;
}

double Planner::getD(int lane) const
{
  double d = (0.5 + lane) * m_lane_width;

  // add a safety margin for the outer lanes
  if (lane == 0) {
    d += m_lane_width * 0.05;
  }
  else if (lane == m_lanes-1) {
    d -= m_lane_width * 0.05;
  }
  return d;
}