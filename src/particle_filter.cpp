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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 1000;
	
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i=0; i<num_particles; i++){
		Particle _particle;
		_particle.id = i;
		_particle.x = dist_x(gen);
		_particle.y = dist_y(gen);
		_particle.theta = dist_theta(gen);	
		_particle.weight = 1.0;
		particles.push_back(_particle);
		weights.push_back(_particle.weight);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen_pred;
	normal_distribution<double> pos_x_pred(0, std_pos[0]);
	normal_distribution<double> pos_y_pred(0, std_pos[1]);
	normal_distribution<double> theta_pred(0, std_pos[2]);
	
	if (fabs(yaw_rate) > 0.0001){
		for (int i=0; i<num_particles; i++){
			particles[i].x = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + pos_x_pred(gen_pred);
			particles[i].y = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + pos_y_pred(gen_pred);
			particles[i].theta = particles[i].theta + yaw_rate*delta_t + theta_pred(gen_pred);
		}
	}
	else{
		for (int i=0; i<num_particles; i++){
			particles[i].x = particles[i].x + velocity*cos(particles[i].theta)*delta_t + pos_x_pred(gen_pred);
			particles[i].y = particles[i].y + velocity*sin(particles[i].theta)*delta_t + pos_y_pred(gen_pred);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i=0; i<observations.size(); i++){
		float dist_min = (observations[i].x-predicted[0].x)*(observations[i].x-predicted[0].x) + (observations[i].y-predicted[0].y)*(observations[i].y-predicted[0].y) ;
		int min_index = 0;
		for (int j=1; j<predicted.size(); j++){
			float dist = (observations[i].x-predicted[j].x)*(observations[i].x-predicted[j].x) + (observations[i].y-predicted[j].y)*(observations[i].y-predicted[j].y) ;
			if (dist < dist_min){
				dist_min = dist;
				min_index = j;
			}
		}
		observations[i].id = min_index;
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
	
	for (int i=0; i<num_particles; i++){
		std::vector<LandmarkObs> observations_in_world;
		for (int j=0; j< observations.size(); j++){
			LandmarkObs _landmarkObs;
			_landmarkObs.x = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
			_landmarkObs.y = particles[i].y + observations[j].y * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);
			_landmarkObs.id = observations[j].id;
			observations_in_world.push_back(_landmarkObs);
		}
		
		std::vector<LandmarkObs> predicted;
		for (int j=0; j< map_landmarks.size(); j++){
			LandmarkObs _landmarkObs;
			_landmarkObs.x = map_landmarks[j].x_f;
			_landmarkObs.y = map_landmarks[j].y_f;
			_landmarkObs.id = map_landmarks[j].id_i;
			predicted.push_back(_landmarkObs);
		}
		dataAssociation(predicted, observations_in_world);
		double prob = 1.0;
		for (int j=0; j< observations.size(); j++){
			double delx = observations[j].x - predicted[observations[j].id].x;
			double dely = observations[j].y - predicted[observations[j].id].y;
			prob *= 1.0/(2*3.14*std_landmark[0]*std_landmark[1])*exp( -( delx*delx/(2*std_landmark[0]*std_landmark[0]) + dely*dely/(2*std_landmark[1]*std_landmark[1]) ) );
		}
		particles[i].weight = prob;
		weights[i] = prob;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	double beta = 0.0;
	int index = rand()%num_particles;
	
	std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);
  
	max_w = weights[0];
	for (int i = 1; i < num_particles; i++){
		if (max_w < weights[i])
			max_w = weights[i];
	}
	std::vector<Particle> particles_temp;
	for (int i = 1; i < num_particles; i++){
		beta += distribution(generator) * 2.0 * max_w;
		while (beta > weights[index]){
			beta -= weights[index]
			index = (index + 1) % num_particles
		}
		particles_temp.push_back(particles[index]);
	}
	particles = particles_temp;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
