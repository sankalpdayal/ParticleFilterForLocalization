/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random> // Need this for sampling from distributions
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.00001 
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	particles.resize(num_particles);
	weights.resize(num_particles);
	
	double init_weight = 1.0;
	
	for (int i=0; i<num_particles; i++){
		particles[i].id = i;
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);	
		particles[i].weight = init_weight;
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	
	normal_distribution<double> pos_x_pred(0.0, std_pos[0]);
	normal_distribution<double> pos_y_pred(0.0, std_pos[1]);
	normal_distribution<double> theta_pred(0.0, std_pos[2]);
	
	for (int i=0; i<num_particles; i++){
		if (fabs(yaw_rate) > 0.0001){
			particles[i].x += velocity/yaw_rate*( sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity/yaw_rate*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta));
			particles[i].theta += yaw_rate*delta_t;
		}
		else{
			particles[i].x += velocity*cos(particles[i].theta)*delta_t;
			particles[i].y += velocity*sin(particles[i].theta)*delta_t;
		}
		particles[i].x += pos_x_pred(gen);
		particles[i].y += pos_y_pred(gen);
		particles[i].theta += theta_pred(gen);	
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (unsigned int i=0; i<observations.size(); i++){
		double dist_min = dist(observations[i].x, observations[i].y, predicted[0].x, predicted[0].y) ;
		int min_index_id = predicted[0].id;
		for (unsigned int j=1; j<predicted.size(); j++){
			double dist_test = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y)  ;
			if (dist_test < dist_min){
				dist_min = dist_test;
				min_index_id = predicted[j].id;
			}
		}
		observations[i].id = min_index_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
	
	double sum_weights = 0.0;
	double sigmaxx = std_landmark[0]*std_landmark[0];
	double sigmayy = std_landmark[1]*std_landmark[1];
		
	for (int i=0; i<num_particles; i++){
		std::vector<LandmarkObs> observations_in_world;
		observations_in_world.resize(observations.size());
		for (unsigned int j=0; j< observations.size(); j++){
			observations_in_world[j].x = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
			observations_in_world[j].y = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);
			observations_in_world[j].id = observations[j].id;
		}
		
		std::vector<LandmarkObs> predicted;
		for (unsigned int j=0; j< map_landmarks.landmark_list.size(); j++){
			double x = map_landmarks.landmark_list[j].x_f;
			double y = map_landmarks.landmark_list[j].y_f;
			int id = map_landmarks.landmark_list[j].id_i;
			double range = dist(x, y, particles[i].x, particles[i].y);
			if (range < sensor_range){
				predicted.push_back(LandmarkObs{ id, x, y });
			}
		}
		dataAssociation(predicted, observations_in_world);
		
		
		vector<int> associations;
		vector<double> sense_x;
		vector<double> sense_y;
		double prob = 1.0;
		
		for (unsigned int j=0; j< observations_in_world.size(); j++){
			double mu_x, mu_y;
			// get the x,y coordinates of the landmark
			for (unsigned int k = 0; k < predicted.size(); k++) {
				if (predicted[k].id == observations_in_world[j].id) {
				  mu_x = predicted[k].x;
				  mu_y = predicted[k].y;
				}
			}
			double delx = observations_in_world[j].x - mu_x;
			double dely = observations_in_world[j].y - mu_y;

			double lm_weight = exp( -( delx*delx/(2*sigmaxx) + (dely*dely/(2*sigmayy)) ) );
			if (lm_weight == 0) {
				prob *= EPS;
			}else {
				prob *= lm_weight;
			}
			associations.push_back(observations_in_world[j].id);
			sense_x.push_back(observations_in_world[j].x);
			sense_y.push_back(observations_in_world[j].y);
		}
		particles[i].weight = prob;
		particles[i] = SetAssociations(particles[i], associations, sense_x, sense_y);
		sum_weights += prob;
	}
	for (int i=0; i<num_particles; i++){
		particles[i].weight /= sum_weights; //coefficient 1.0/(2*MI_PI*std_landmark[0]*std_landmark[1]) will cancel here
		weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	discrete_distribution<int> dist_particles(weights.begin(), weights.end());
	std::vector<Particle> particles_temp;
	particles_temp.resize(num_particles);
	for (int i = 0; i < num_particles; i++){
		auto index = dist_particles(gen);
		particles_temp[i]= std::move(particles[index]);
	}
	particles = std::move(particles_temp);
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
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
