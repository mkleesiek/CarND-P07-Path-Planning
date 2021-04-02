#ifndef PATH_PLANNING_PLANNER_H
#define PATH_PLANNING_PLANNER_H

#include <vector>
#include <deque>
#include <cstdint>

/**
 * A path planner for the Udacity self-driving car path planning project.
 *
 * To be used in conjunction with the Term 3 Simulator.
 * @see https://github.com/udacity/self-driving-car-sim/releases/tag/T3_v1.2
 */
class Planner {
public:
  /// Holding coordinates of the simulator track, allowing conversion between cartesian and frenet coordinates.
  struct Map {
    std::vector<double> waypoints_x;
    std::vector<double> waypoints_y;
    std::vector<double> waypoints_s;
    double max_s;
  };

  /// Holding the state of a vehicle (s and d coordinate, their velocities and the time).
  struct CarState {
    double s;
    double d;
    double s_dot;
    double d_dot;
    double t;
  };

  /// A cartesian coordinate including time.
  struct PathPoint {
    double x;
    double y;
    double t;
  };

  /**
   * Reset the planner to it's initialization state.
   */
  void reset();

  /**
   * Initialize map data from the simulator.
   * This function might interpolate the input data to create a more fine-grained and smooth internal representation.
   */
  void setMapData(const std::vector<double>& map_x, const std::vector<double>& map_y, const std::vector<double>& map_s,
    double max_s);

  /**
   * Set the speed limit for path planning.
   *
   * @param mph speed limit in miles-per-hours.
   */
  void setSpeedLimitMPH(double mph);

  /**
   * Update the own car coordinates from simulator input data.
   * The vehicle's time relative to the already planned path is resolved based on it's current position.
   */
  void updateCarLocalization(double x, double y, double s, double d, double yaw, double speedMPH);

  /**
   * Update the internal list of other vehicles on the track, based on sensor fusion data provided by the simulator.
   */
  void updateSensorFusion(const std::vector<std::vector<double>>& sensor_data);

  /**
   * Plan the upcoming vehicle trajectory and provide it in arrays of cartesian x/y coordinates.
   * @param[output] next_x
   * @param[output] next_y
   */
  void planTrajectory(std::vector<double>& next_x, std::vector<double>& next_y);

private:
  /// Suggest the next waypoint
  CarState suggestNextWaypoint(const CarState& previous) const;
  /// Generate a jerk-free trajectory using Quintic Polynomial solving
  std::vector<CarState> generateJerkMinimizingTrajectory(const std::deque<CarState>& waypoints) const;
  /// Check for collisions with other vehicles based on our current and suggested next waypoint
  std::pair<bool, double> checkCollision(const CarState& me_current, const CarState& me_next) const;

  /// Convert to Frenet coordinates
  std::pair<double, double> getFrenet(double x, double y, double theta) const;
  /// Convert to Cartesion coordinates
  std::pair<double, double> getXY(double s, double d) const;

  /// Normalize the Frenet s coordinate based on our map's maximum s
  double normalize_s(double s) const;
  /// Determine the distance between two s positions on our closed track
  double distance_s(double s_front, double s_back) const;

  /// Get lane number for a given d coordinate
  int getLane(double d) const;
  /// Get the d position of a lane center
  double getD(int lane) const;

  /// Speed limit in m/s
  double m_speedLimit = 10.0;

  /// Number of lanes per direction
  int m_lanes = 3;
  /// Lane width in m
  double m_lane_width = 4.0;

  /// update interval in s
  double m_update_interval = 0.02;

  /// duration for looking (planing ahead) in s
  double m_look_ahead = 4.0;
  /// minimum duration of planned path (in s) which will not be discarded or updated anymore
  double m_buffer = 0.5;
  /// minimum distance to other vehicels in m
  double m_distance_to_others = 30.0;

  /// Map data
  Map m_map{};
  /// Our car's current state
  CarState m_current_state{};
  /// States of all other vehicles
  std::vector<CarState> m_other_cars;
  /// Currently planned waypoints in Frenet coordinates (each a couple of seconds apart)
  std::deque<CarState> m_planned_waypoints;
  /// Currently planned path in Cartesian coordinates
  std::vector<PathPoint> m_planned_path;
};

#endif //PATH_PLANNING_PLANNER_H
