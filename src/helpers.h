#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>
#include <string>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Dense>

// for convenience
using std::string;
using std::vector;
using Eigen::VectorXd;
using Eigen::MatrixXd;

/**
 * Conversion from miles-per-hour to meter-per-second.
 *
 * @param mph speed in miles-per-hour.
 * @return speed in meters-per-second
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
 * @return an array of length 6, each value corresponding to a coefficent in
 *   the polynomial:
 *   s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5
 */
inline vector<double> jmt(const vector<double>& start, const vector<double>& end, double T)
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

/**
 * Check if two floating point numbers are equal within a given precision.
 * @tparam T Floating point type.
 * @param a 1st number.
 * @param b 2nd number.
 * @param epsilon Precision.
 * @return True if 'a' and 'b' are within epsilon.
 */
template<class T>
constexpr bool essentially_equal(T a, T b, T epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

/**
 * Calculate a polynomial.
 * @param x Parameter value.
 * @param coeff Array of polynomial coefficients.
 * @return
 */
inline double polynom(double x, const vector<double>& coeff)
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

/**
 * Calculate the first derivative of a polynomial.
 * @param x Parameter value.
 * @param coeff Array of polynomial coefficients.
 * @return
 */
inline double polynom_gradient(double x, const vector<double>& coeff)
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

//
// Helper functions provided by Udacity down below:
//

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
//   else the empty string "" will be returned.
inline string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

//
// Helper functions related to waypoints and converting from XY to Frenet
//   or vice versa
//

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
constexpr double deg2rad(double x) { return x * pi() / 180; }
constexpr double rad2deg(double x) { return x * 180 / pi(); }

// Calculate squared distance between two points
inline double distance_squared(double x1, double y1, double x2, double y2) {
  const double x_diff = (x2-x1);
  const double y_diff = (y2-y1);
  return x_diff*x_diff+y_diff*y_diff;
}

// Calculate distance between two points
inline double distance(double x1, double y1, double x2, double y2) {
  return sqrt(distance_squared(x1, y1, x2, y2));
}

// Calculate closest waypoint to current x, y position
inline int ClosestWaypoint(double x, double y, const vector<double> &maps_x,
                    const vector<double> &maps_y) {
  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); ++i) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x,y,map_x,map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

// Returns next waypoint of the closest waypoint
inline int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x,
                 const vector<double> &maps_y) {
  int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y-y),(map_x-x));

  double angle = fabs(theta-heading);
  angle = std::min(2*pi() - angle, angle);

  if (angle > pi()/2) {
    ++closestWaypoint;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
inline vector<double> getFrenet(double x, double y, double theta,
                         const vector<double> &maps_x, 
                         const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

  int prev_wp;
  prev_wp = next_wp-1;
  if (next_wp == 0) {
    prev_wp  = maps_x.size()-1;
  }

  double n_x = maps_x[next_wp]-maps_x[prev_wp];
  double n_y = maps_y[next_wp]-maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
  double proj_x = proj_norm*n_x;
  double proj_y = proj_norm*n_y;

  double frenet_d = distance(x_x,x_y,proj_x,proj_y);

  //see if d value is positive or negative by comparing it to a center point
  double center_x = 1000-maps_x[prev_wp];
  double center_y = 2000-maps_y[prev_wp];
  double centerToPos = distance(center_x,center_y,x_x,x_y);
  double centerToRef = distance(center_x,center_y,proj_x,proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; ++i) {
    frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
  }

  frenet_s += distance(0,0,proj_x,proj_y);

  return {frenet_s,frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
inline vector<double> getXY(double s, double d, const vector<double> &maps_s,
                     const vector<double> &maps_x, 
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1))) {
    ++prev_wp;
  }

  int wp2 = (prev_wp+1)%maps_x.size();

  double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),
                         (maps_x[wp2]-maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s-maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
  double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

  double perp_heading = heading-pi()/2;

  double x = seg_x + d*cos(perp_heading);
  double y = seg_y + d*sin(perp_heading);

  return {x,y};
}

#endif  // HELPERS_H