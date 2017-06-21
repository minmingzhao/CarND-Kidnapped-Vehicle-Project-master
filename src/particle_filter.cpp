/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>

#include "particle_filter.h"

using namespace std;
default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;


	normal_distribution<double> dist_x_noise_init(0, std[0]);
	normal_distribution<double> dist_y_noise_init(0, std[1]);
	normal_distribution<double> dist_theta_noise_init(0, std[2]);

	//init particles
	for (int i = 0; i<num_particles; i++) {
		Particle p;
		p.id = i;
		p.x = x;
		p.y = y;
		p.theta = theta;
		p.weight = 1.0;

		//add noise term
		p.x += dist_x_noise_init(gen);
		p.y += dist_y_noise_init(gen);
		p.theta += dist_theta_noise_init(gen);

		particles.push_back(p);
	}

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double>  dist_x_noise(0, std_pos[0]);
	normal_distribution<double>  dist_y_noise(0, std_pos[1]);
	normal_distribution<double>  dist_theta_noise(0, std_pos[2]);

	for (int i = 0; i < num_particles; i++) {
		if (fabs(yaw_rate) < 0.00001) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta)+dist_x_noise(gen);
			particles[i].y += velocity * delta_t * sin(particles[i].theta)+dist_y_noise(gen);
		}
		else {
			particles[i].x += velocity/yaw_rate * (sin(particles[i].theta+yaw_rate*delta_t) - sin(particles[i].theta))+dist_x_noise(gen);
			particles[i].y += velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t))+dist_y_noise(gen);

		}
		particles[i].theta += yaw_rate * delta_t + dist_theta_noise(gen);

	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (unsigned int i =0; i < observations.size(); i++) {
		LandmarkObs obs = observations[i];

		double min_dist = numeric_limits<double>::max();

		int map_id = -1;

		for (unsigned int j=0; j<predicted.size();j++) {
			LandmarkObs p = predicted[j];

			double cur_dist = dist(obs.x, obs.y, p.x, p.y);
			if (cur_dist < min_dist) {
				min_dist = cur_dist;
				map_id = p.id;
			}

		}
		//associate this observation to that map_id landmark
		observations[i].id = map_id;
		//cout<<"Landmark Index: " <<map_id<< "is found for observation (x,y)=("<< obs.x <<", "<<obs.y <<")"<<endl;
	}
}


double prob(double x, double y, float xm, float ym, double std_landmark[]) {
	double xpart = -0.5*(x - xm)*(x - xm) / (std_landmark[0] * std_landmark[0]);
	double ypart = -0.5*(y - ym)*(y - ym) / (std_landmark[1] * std_landmark[1]);
	return exp(xpart + ypart) / (2 * M_PI*std_landmark[0] * std_landmark[1]);
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	//Map temp_landmarks;
	//temp_landmarks.landmark_list.clear(); // clear the vector
	weights.clear(); // clear the previous weights vector

	for (int i = 0; i < num_particles; i++)
	{
		double ctheta = cos(particles[i].theta);
		double stheta = sin(particles[i].theta);
		double p_x = particles[i].x;
		double p_y = particles[i].y;

		//temp_landmarks.landmark_list.clear(); // clear the vector
		vector<LandmarkObs> predictions;
		predictions.clear();
		for (int z = 0; z < map_landmarks.landmark_list.size(); ++z)
		{
            if (fabs(p_x-map_landmarks.landmark_list[z].x_f)<sensor_range && fabs(p_y-map_landmarks.landmark_list[z].y_f)<sensor_range) {
            //if (dist(p_x, p_y, map_landmarks.landmark_list[z].x_f, map_landmarks.landmark_list[z].y_f) < sensor_range) {
            	predictions.push_back(LandmarkObs{map_landmarks.landmark_list[z].id_i, map_landmarks.landmark_list[z].x_f, map_landmarks.landmark_list[z].y_f});
			}
		}

		particles[i].weight = 1.0; // set it to 1 before the next calculations begin

		//LandmarkObs LMObs;
		//vector<LandmarkObs> transform_obs;
		//transform_obs.clear();
		for (unsigned int j=0; j < observations.size(); j++)
		{
			double t_x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + p_x;
			double t_y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + p_y;
			//transform_obs.push_back(LandmarkObs{observations[j].id, temp_x, temp_y});
			//cout << particle.id << "local" << LMObs.x << " " << LMObs.y << "global" << temp_x << " " << temp_y;
			//particle.sense_x[j] = temp_x; //obs rotated then translated
			//particle.sense_y[j] = temp_y; //obs rotated then translated

			double min = 999.9;
			int min_id=0;
			for (int k = 0; k < predictions.size(); ++k)
			{
				double distance = dist(t_x, t_y, predictions[k].x, predictions[k].y);
				if (distance < min)
				{
					min = distance;
					min_id = predictions[k].id;
				}
			}

			if (observations.size() == 0) {
				particles[i].weight = 0.0;
			}
			else
			{
				particles[i].weight *= prob(t_x, t_y, map_landmarks.landmark_list[min_id-1].x_f, map_landmarks.landmark_list[min_id-1].y_f, std_landmark);
			}
			//cout << " weight: " << particle.weight << "\n";

		}
		weights.push_back(particles[i].weight);
		predictions.clear();
	}
}
void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particles;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d(weights.begin(), weights.end());
	std::map<int, int> m;
	for (int n = 0; n<num_particles; ++n) {
		new_particles.push_back(particles[d(gen)]);
	}
	particles = new_particles; // resampled set of particles
	weights.clear(); // clear the previous weights vector
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
