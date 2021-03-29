#ifndef PATH_PLANNING_PLANNER_H
#define PATH_PLANNING_PLANNER_H

#include <vector>
#include <deque>
#include <unordered_map>

class Planner {
public:
  struct Map {
    std::vector<double> waypoints_x;
    std::vector<double> waypoints_y;
    std::vector<double> waypoints_s;
    std::vector<double> waypoints_dx;
    std::vector<double> waypoints_dy;
    double max_s;
  };

  struct CarState {
    double x;
    double y;
    double s;
    double d;
    double s_dot;
    double d_dot;
  };

  void setMapData(const Map& map);
  void setSpeedLimitMPH(double mph);

  void updateCarLocalization(double x, double y, double s, double d, double yaw, double speedMPH);
  void updateSensorFusion(const std::vector<std::vector<double>>& sensor_data);
  void updatePreviousPath(const std::vector<double>& path_x, const std::vector<double>& path_y);

  CarState suggestTrajectory(const CarState& previous);
  void suggestPath(std::vector<double>& next_x, std::vector<double>& next_y);

private:
  std::pair<double, double> getFrenet(double x, double y, double theta) const;
  std::pair<double, double> getXY(double s, double d) const;

  /// Map data
  Map m_map;
  /// Speed limit in m/s
  double m_speedLimit;

  int m_lanes = 3;
  double m_lane_width = 4.0;

  double m_update_interval = 0.02;
  double m_max_accel = 10.0;

  double m_time_per_waypoint = 3.0;
  double m_buffer = 1.0;

  CarState m_current_state;
  std::deque<CarState> m_previous_states;
  std::unordered_map<int, CarState> m_other_cars;
};

#endif //PATH_PLANNING_PLANNER_H
