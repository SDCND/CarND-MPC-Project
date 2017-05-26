#define _USE_MATH_DEFINES
#include <cmath>

#include "uWS/uWS.h"
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "MPC.h"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"


// for convenience
using json = nlohmann::json;
using namespace std;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::stringstream hasData(std::string s) 
{
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_last_of("]");

  if (found_null != std::string::npos) 
	{
    return stringstream();
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) 
	{
	  std::stringstream tmp = std::stringstream();
	  tmp.str(s.substr(b1, b2 - b1 + 1));
	  return tmp;
  }

  return std::stringstream();
}


// polinomial evaluation
double polyeval(Eigen::VectorXd coeffs, double x) 
{
		double result = 0.0;

		for (int i = 0; i < coeffs.size(); i++) 
		{
				result += coeffs[i] * pow(x, i);
		}

		return result;
}


// polynomial fitting, adapted from https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,	int order) 
{
		assert(xvals.size() == yvals.size());
		assert(order >= 1 && order <= xvals.size() - 1);

		Eigen::MatrixXd A(xvals.size(), order + 1);

		for (int i = 0; i < xvals.size(); i++) 
		{
				A(i, 0) = 1.0;
		}

		for (int j = 0; j < xvals.size(); j++) 
		{
				for (int i = 0; i < order; i++) 
				{
						A(j, i + 1) = A(j, i) * xvals(j);
				}
		}

		auto Q = A.householderQr();
		auto result = Q.solve(yvals);

		return result;
}


int main(int argc, const char *argv[])
{
  uWS::Hub h;
	
	MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) 
	{
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {
	  auto s = hasData(std::string(data));

      if (s.str() != "") 
			{
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();

        if (event == "telemetry") 
				{
						// j[1] is the data JSON object
						vector<double> ptsx = j[1]["ptsx"];
						vector<double> ptsy = j[1]["ptsy"];

						double px = j[1]["x"];
						double py = j[1]["y"];
						double psi = j[1]["psi"];
						double v = j[1]["speed"];

						/*
						* TODO: Calculate steeering angle and throttle using MPC.
						*
						* Both are in between [-1, 1].
						*
						*/

						vector<double> next_x_vals;
						vector<double> next_y_vals;

						Eigen::VectorXd ptsx_X(ptsx.size());
						Eigen::VectorXd ptsy_X(ptsy.size());

						for (int i = 0; i < ptsx.size(); i++)
						{
								double diff_x = ptsx[i] - px;
								double diff_y = ptsy[i] - py;

								// transformation
								ptsx_X(i) = diff_x * cos(psi) + diff_y * sin(psi);
								ptsy_X(i) = -diff_x * sin(psi) + diff_y * cos(psi);

								next_x_vals.push_back(ptsx_X(i));
								next_y_vals.push_back(ptsy_X(i));
						}

						auto coeffs = polyfit(ptsx_X, ptsy_X, 3); 

						// cross track error
						double cte = polyeval(coeffs, 0.0);

						// orientation error
						double epsilon = atan(coeffs[1]);

						Eigen::VectorXd state(6);
						state << 0.0, 0.0, 0.0, v, cte, epsilon;
						std::cout << "cte: " << cte	<< ", psi: " << psi << ", epsilon: " << epsilon << std::endl;
						
						auto vars = mpc.Solve(state, coeffs);
						double steer_value = -vars[0];
						double throttle_value = vars[1];

						json msgJson;
						msgJson["steering_angle"] = steer_value;
						msgJson["throttle"] = throttle_value;


						//Display the MPC predicted trajectory 
						vector<double> mpc_x_vals;
						vector<double> mpc_y_vals;

						//.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
						// the points in the simulator are connected by a Green line

						msgJson["mpc_x"] = mpc_x_vals;
						msgJson["mpc_y"] = mpc_y_vals;

						//.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
						// the points in the simulator are connected by a Yellow line

						msgJson["next_x"] = next_x_vals;
						msgJson["next_y"] = next_y_vals;

						auto msg = "42[\"steer\"," + msgJson.dump() + "]";
						std::cout << msg << std::endl;

						// Latency
						// The purpose is to mimic real driving conditions where
						// the car does actuate the commands instantly.
						//
						// Feel free to play around with this value but should be to drive
						// around the track with 100ms latency.
						//
						// NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
						// SUBMITTING.
						this_thread::sleep_for(chrono::milliseconds(100));
						ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } 
			else 
			{
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) 
	{
    const std::string s = "<h1>Hello world!</h1>";

    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) 
	{
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) 
	{
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;

  if (h.listen("0.0.0.0", port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }

  h.run();
}
