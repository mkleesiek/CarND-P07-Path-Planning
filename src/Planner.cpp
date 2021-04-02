#include "Planner.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "helpers.h"
#include "spline.h"
#include <iostream>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace {
  /**
   * Conversion from miles-per-hour to meter-per-second.
   */
  constexpr double mph_to_ms(double mph)
  {
    return mph * 1609.34 / 3600.0;
  }

  /**
   * Calculate the Jerk Minimizing Trajectory that connects the initial state
   * to the final state in time T.
   *
   * @param start - the vehicles start location given as a length three array
   *   corresponding to initial values of [s, s_dot, s_double_dot]
   * @param end - the desired end state for vehicle. Like "start" this is a
   *   length three array.
   * @param T - The duration, in seconds, over which this maneuver should occur.
   *
   * @output an array of length 6, each value corresponding to a coefficent in
   *   the polynomial:
   *   s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
   */
  vector<double> jmt(const vector<double>& start, const vector<double>& end, double T)
  {
    VectorXd alpha_012(3);
    alpha_012 << start[0], start[1], 0.5*start[2];

    VectorXd s_diff(3);
    s_diff << end[0] - (start[0] + start[1]*T + 0.5*start[2]*T*T),
              end[1] - (start[1] + start[2]*T),
              end[2] - start[2];

    MatrixXd matrix_t(3, 3);
    matrix_t << T*T*T,   T*T*T*T,   T*T*T*T*T,
                3.0*T*T, 4.0*T*T*T, 5.0*T*T*T*T,
                6.0*T,   12.0*T*T,  20.0*T*T*T;

    const VectorXd alpha_345 = matrix_t.inverse() * s_diff;

    return {alpha_012(0), alpha_012(1), alpha_012(2), alpha_345(0), alpha_345(1), alpha_345(2)};
  }

  double polynom(double x, const vector<double>& coeff)
  {
    if (coeff.size() < 1) {
      return 0.0;
    }

    double y = coeff[0];
    double x_p = 1.0;
    for (size_t i = 1; i < coeff.size(); i++) {
      x_p *= x;
      y += coeff[i] * x_p;
    }

    return y;
  }

  double polynom_gradient(double x, const vector<double>& coeff)
  {
    if (coeff.size() < 2) {
      return 0.0;
    }

    double y = coeff[1];
    double x_p = 1.0;
    for (size_t i = 2; i < coeff.size(); i++) {
      x_p *= x;
      y += coeff[i] * i * x_p;
    }

    return y;
  }

  template<class T>
  constexpr bool essentially_equal(T a, T b, T epsilon)
  {
      return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
  }
}

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
  for (size_t k = 0; k < 5; k++) {
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
  // condition for a clean start
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
//    cout << "adding new waypoint" << endl;
  }

  // generate an interpolated jerk-free trajectory from the suggested waypoints
  auto planned_jerkfree = generateJerkfreeTrajectory(m_planned_waypoints);

//  // normalize s coordinates
//  for (auto& waypoint : m_planned_waypoints) {
//    waypoint.s = normalize_s(waypoint.s);
//  }

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
    if (t < m_current_state.t + m_update_interval/2.0) {
      continue;
    }

    PathPoint p{sx(t), sy(t), t};

    next_x.push_back(p.x);
    next_y.push_back(p.y);
    m_planned_path.push_back(std::move(p));
  }
}

vector<Planner::CarState> Planner::generateJerkfreeTrajectory(const deque<CarState>& waypoints) const
{
  vector<CarState> jerkfree_traj;
  jerkfree_traj.push_back(waypoints.front());

  for (size_t i = 1; i < waypoints.size(); i++) {

    const auto& last_state = waypoints[i-1];
    auto next_state = waypoints[i];

    const double traj_duration = next_state.t - last_state.t;
    const auto coeff_s = jmt({last_state.s, last_state.s_dot, 0.0}, {next_state.s, next_state.s_dot, 0.0}, traj_duration);
    const auto coeff_d = jmt({last_state.d, last_state.d_dot, 0.0}, {next_state.d, next_state.d_dot, 0.0}, traj_duration);

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
      jerkfree_traj.push_back(inter_state);
    }
  }

  return jerkfree_traj;
}

Planner::CarState Planner::suggestNextWaypoint(const CarState& previous) const
{
  // segment duration when staying in same lane
  constexpr double dt_straight = 2.0;
  // segment duration for a lane change
  constexpr double dt_lanechange = 3.0;

  const int current_lane = getLane(previous.d);

  CarState next_straight{};
  next_straight.t = previous.t + dt_straight;
  next_straight.s_dot = min(previous.s_dot + m_max_accel * dt_straight, m_speedLimit);
  next_straight.s = previous.s + (next_straight.s_dot + previous.s_dot) / 2.0 * dt_straight;
  next_straight.d = getD(current_lane);
  auto collision_straight = checkCollision(previous, next_straight);

  CarState next_left{};
  next_left.t += previous.t + dt_lanechange;
  next_left.s_dot = min(previous.s_dot + m_max_accel / 2.0 * dt_lanechange, m_speedLimit);
  next_left.s = previous.s + (next_left.s_dot + previous.s_dot) / 2.0 * dt_lanechange;
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

  // 2) keep current lane otherwise,
  //    or if we have just gotten off the start
  if (!collision_straight.first || previous.t < 5.0) {
    return next_straight;
  }

  // 3) we'd have to slow down in order not to collide, so consider changing lanes

  // check left-hand lane
  if (current_lane > 0 && !collision_left.first) {
    return next_left;
  }
  // check right-hand lane
  else if (current_lane < m_lanes-1 && !collision_right.first) {
    return next_right;
  }


  // 4) we can't change langes right now and have to reduce speed to avoid collision

  // adjust to the spped of the car in front
  next_straight.s_dot = min(next_straight.s_dot, collision_straight.second * 0.9);
  // but don't exceed the (negative) acceleration limit
  next_straight.s_dot = max(previous.s_dot - m_max_accel * dt_straight, next_straight.s_dot);
  // estimate distance based on target velocity
  next_straight.s = previous.s + (next_straight.s_dot + previous.s_dot) / 2.0 * dt_straight;

  return next_straight;
}

pair<bool, double> Planner::checkCollision(const CarState& me_current, const CarState& me_next) const
{
  const int my_lane_now = getLane(me_current.d);
  const int my_lane_next = getLane(me_next.d);

  double min_distance = m_distance_to_others;
  CarState candidate{};

  for (const auto& other : m_other_cars) {

    const int others_lane = getLane(other.d);

    // skip vehicles which will not be on my lane
    if (my_lane_next != others_lane) {
      continue;
    }

    // skip other vehicles which are right behind me on my lane!
    if (others_lane == my_lane_now && me_current.t >= other.t && normalize_s(me_current.s) > other.s) {
      continue;
    }

    // extrapolate the other car position to my next one
    double s_other = other.s + other.s_dot * (me_next.t - other.t);
    double distance_next = s_other - normalize_s(me_next.s);

    // check for collision on my next state
    if (fabs(distance_next) < min_distance) {
      candidate = other;
      min_distance = fabs(distance_next);
    }

    // extrapolate the other car position to my current one, if on the same lane and ahead
    if (my_lane_now == my_lane_next) {
      s_other = other.s + other.s_dot * (me_current.t - other.t);
      const double distance_current = s_other - normalize_s(me_current.s);

      // check for collision on my current state
      if (distance_current > 0.0 && distance_current < min_distance) {
        candidate = other;
        min_distance = distance_current;
      }
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
  s = fmod(s, m_map.max_s);
  const auto xy = ::getXY(s, d, m_map.waypoints_s, m_map.waypoints_x, m_map.waypoints_y);
  return make_pair(xy[0], xy[1]);
}

double Planner::normalize_s(double s) const
{
  return fmod(s, m_map.max_s);
}

int Planner::getLane(double d) const
{
  return d / m_lane_width;
}

double Planner::getD(int lane) const
{
  double d = (0.5 + lane) * m_lane_width;
  if (lane == 0) {
    d += m_lane_width * 0.05;
  }
  else if (lane == m_lanes-1) {
    d -= m_lane_width * 0.05;
  }
  return d;
}