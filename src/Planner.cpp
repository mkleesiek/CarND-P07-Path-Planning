#include "Planner.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "helpers.h"
#include "spline.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace {
  constexpr double mph_to_ms(double mph)
  {
    return mph * 1609.34 / 3600.0;
  }

  vector<double> jmt(const vector<double>& start, const vector<double>& end, double T)
  {
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
   *
   * EXAMPLE
   *   > JMT([0, 10, 0], [10, 10, 0], 1)
   *     [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
   */

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
    if (coeff.empty()) {
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
}

void Planner::setMapData(const Map& map)
{
  m_map = map;
}

void Planner::setSpeedLimitMPH(double mph)
{
  m_speedLimit = mph_to_ms(mph);
}

void Planner::updateCarLocalization(double x, double y, double s, double d, double yaw, double speedMPH)
{
  const double speed = mph_to_ms(speedMPH);
  const double x_dot = speed * cos(deg2rad(yaw));
  const double y_dot = speed * sin(deg2rad(yaw));

  const auto sd_2 = getFrenet(x + x_dot * 1.0, y + y_dot * 1.0, deg2rad(yaw));
  const double s_dot = sd_2.first - s / 1.0;
  const double d_dot = sd_2.second - d / 1.0;

  m_current_state = {x, y, s, d, s_dot, d_dot};
}

void Planner::updateSensorFusion(const vector<vector<double>>& sensor_data)
{
  m_other_cars.clear();

  for (const auto& dataSet : sensor_data)
  {
    const int id = (int) dataSet[0];
    const double s = dataSet[5];
    const double d = dataSet[6];
    const double heading = atan2(dataSet[4], dataSet[3]);
    const auto sd_2 = getFrenet(dataSet[1] + dataSet[3] * 1.0, dataSet[2] + dataSet[4] * 1.0, heading);
    const double s_dot = sd_2.first - s / 1.0;
    const double d_dot = sd_2.second - d / 1.0;

    m_other_cars.emplace(id, CarState{dataSet[1], dataSet[2], s, d, s_dot, d_dot});
  }
}

void Planner::updatePreviousPath(const std::vector<double>& path_x, const std::vector<double>& path_y)
{
  m_previous_states.clear();

  const size_t n_previous = min(path_x.size(), path_y.size());
  for (size_t i = 0; i < n_previous; i++) {
    m_previous_states.push_back(CarState{path_x[i], path_y[i]});
  }
}

void Planner::suggestPath(vector<double>& next_x, vector<double>& next_y)
{
  vector<double> t_vals, x_vals, y_vals;
  t_vals.push_back(0.0);
  x_vals.push_back(m_current_state.x);
  y_vals.push_back(m_current_state.y);

  const size_t n_prev = m_previous_states.size();

  double t_ahead = 0.0;
  for (size_t i = 0; i < n_prev; i++) {
    t_ahead = m_update_interval * (i + 1);
    t_vals.push_back(t_ahead);
    x_vals.push_back(m_previous_states[i].x);
    y_vals.push_back(m_previous_states[i].y);

    next_x.push_back(m_previous_states[i].x);
    next_y.push_back(m_previous_states[i].y);
  }

  CarState end = m_current_state;

  if (m_previous_states.size() > 1) {
    auto prev1_state = m_previous_states[n_prev-1];
    auto prev2_state = m_previous_states[n_prev-2];

    double heading = atan2(prev1_state.y-prev2_state.y, prev1_state.x-prev2_state.x);

    auto prev1_sd = getFrenet(prev1_state.x, prev1_state.y, heading);
    auto prev2_sd = getFrenet(prev2_state.x, prev2_state.y, heading);

    double s_dot = (prev1_sd.first - prev2_sd.first) / m_update_interval;
    double d_dot = (prev1_sd.second - prev2_sd.second) / m_update_interval;

    end = CarState{prev1_state.x, prev1_state.y, prev1_sd.first, prev1_sd.second, s_dot, d_dot};
  }

  while (t_ahead < m_buffer || t_vals.size() < 4) {

    CarState next = suggestTrajectory(end);

    const auto coeff_s = jmt({end.s, end.s_dot, 0.0}, {next.s, next.s_dot, 0.0}, m_time_per_waypoint);
    const auto coeff_d = jmt({end.d, 0.0, 0.0}, {next.d, 0.0, 0.0}, m_time_per_waypoint);

    constexpr size_t jmt_iterations = 4;

    for (size_t k = 0; k < jmt_iterations; k++) {
      double dt = m_time_per_waypoint * (k+1) / jmt_iterations;
      next.s = polynom(dt, coeff_s);
      next.d = polynom(dt, coeff_d);

      auto xy = getXY(next.s, next.d);
      next.x = xy.first;
      next.y = xy.second;

      t_vals.push_back(t_ahead + dt);
      x_vals.push_back(next.x);
      y_vals.push_back(next.y);
    }

    t_ahead += m_time_per_waypoint;
    end = next;
  }

  tk::spline sx(t_vals, x_vals);
  tk::spline sy(t_vals, y_vals);

  size_t n_interpolated = t_ahead / m_update_interval;

  for (size_t i = next_x.size(); i < n_interpolated; i++) {
    double t = t_ahead * (i+1) / n_interpolated;
    next_x.push_back(sx(t));
    next_y.push_back(sy(t));
  }
}

Planner::CarState Planner::suggestTrajectory(const CarState& previous)
{
  CarState next{};
  next.s_dot = m_speedLimit;
  next.s_dot = min(previous.s_dot + m_max_accel / 2.0 * m_time_per_waypoint, next.s_dot);

  next.s = previous.s + (next.s_dot + previous.s_dot) / 2.0 * m_time_per_waypoint;
  next.d = (previous.s > 300) ? m_lane_width * (0.5 + 2) : m_lane_width * (0.5 + 1);

  return next;
}

pair<double, double> Planner::getFrenet(double x, double y, double theta) const
{
  const auto sd = ::getFrenet(x, y, theta, m_map.waypoints_x, m_map.waypoints_y);
  return make_pair(sd[0], sd[1]);
}

pair<double, double> Planner::getXY(double s, double d) const
{
  const auto xy = ::getXY(s, d, m_map.waypoints_s, m_map.waypoints_x, m_map.waypoints_y);
  return make_pair(xy[0], xy[1]);
}
