#ifndef PATH_PLANNING_PLANNER_H
#define PATH_PLANNING_PLANNER_H

#include <vector>
#include <deque>
#include <cstdint>

class Planner {
public:
  struct Map {
    std::vector<double> waypoints_x;
    std::vector<double> waypoints_y;
    std::vector<double> waypoints_s;
    double max_s;
  };

  struct CarState {
    double s;
    double d;
    double s_dot;
    double d_dot;
    double t;
  };

  struct PathPoint {
    double x;
    double y;
    double t;
  };

  void setMapData(const std::vector<double>& map_x, const std::vector<double>& map_y, const std::vector<double>& map_s,
    double max_s);
  void setSpeedLimitMPH(double mph);

  void updateCarLocalization(double x, double y, double s, double d, double yaw, double speedMPH);
  void updateSensorFusion(const std::vector<std::vector<double>>& sensor_data);

  void planTrajectory(std::vector<double>& next_x, std::vector<double>& next_y);

private:
  CarState suggestNextWaypoint(const CarState& previous) const;
  std::pair<bool, CarState> checkCollision(const CarState& my_car) const;

  std::pair<double, double> getFrenet(double x, double y, double theta) const;
  std::pair<double, double> getXY(double s, double d) const;

  int getLane(double d) const;
  double getD(int lane) const;

  /// Map data
  Map m_map;
  /// Speed limit in m/s
  double m_speedLimit;

  int m_lanes = 3;
  double m_lane_width = 4.0;

  double m_update_interval = 0.02;
  double m_max_accel = 3.0;

//  double m_trajectory_interval = 2.5;
  double m_buffer = 0.5;

  double m_distance_to_others = 35.0;

  CarState m_current_state = CarState{};
  std::vector<CarState> m_other_cars;
  std::deque<CarState> m_planned_waypoints;
  std::vector<PathPoint> m_planned_path;
};

#endif //PATH_PLANNING_PLANNER_H
