#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "flock.h"
#include "boid.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Flock::Flock(int num_boids) {
    this->num_boids = num_boids;
    buildGrid();
}

Flock::~Flock() {
    boids.clear();
}

void Flock::buildGrid() {
    // populates flock with boids
    
    double numBigBirds = 0.02;
    for (int i = 0; i < num_boids; i++) {
        Vector3D position = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        Vector3D velocity = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        Vector3D acceleration = Vector3D(0, 0, 0);
        bool is_pred = rand() % 100 < numBigBirds * 100;
        Boid boid = Boid(position, velocity, acceleration, is_pred);
        boids.emplace_back(boid);
    }
}

void Flock::simulate(double frames_per_sec, double simulation_steps, FlockParameters *fp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
    double radius;
    double separationRadius;
    double boundaryFactor;
    double separationFactor;
    double alignmentFactor;
    double cohesionFactor;
    double predRange;
    double maxSpeed;
    double minSpeed;
    double bias;
    double hungies; 
    double predTurnFactor;
    bool cohesion;
    bool alignment;
    bool separation;
    

    radius = fp->radius;
    separationRadius = fp->separationRadius;
    boundaryFactor = fp->boundaryFactor;
    separationFactor = fp->separationFactor;
    alignmentFactor = fp->alignmentFactor;
    cohesionFactor = fp->cohesionFactor;
    predRange = fp->predRange;
    hungies = fp ->hungies;
    predTurnFactor = fp->predTurnFactor;

    maxSpeed = fp->maxSpeed;
    minSpeed = fp->minSpeed;
    bias = fp->bias;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;
    
    if (!fp->enable_cohesion) {
        cohesionFactor = 0;
    }
    if (!fp->enable_alignment) {
        alignmentFactor = 0;
    }
    if (!fp->enable_separation) {
        separationFactor = 0;
    }
    
    vector<Boid> left;
    vector<Boid> right;

    for (int i = 0; i < boids.size(); i++) {
        if (i < boids.size() / 2) {
            left.emplace_back(boids[i]);
        } else {
            right.emplace_back(boids[i]);
        }
    }
    
    // iterate through all boids
    
    for (Boid &boid : boids) {
        int neighbors = 0;
        Vector3D averagePosition = Vector3D(0, 0, 0);
        Vector3D averageVelocity = Vector3D(0, 0, 0);
        Vector3D separationVelocity = Vector3D(0, 0, 0);
        
        // iterate through all neighboring boids
        for (Boid &neighbor : boids) {
            if(boid.isPredator || neighbor.isPredator)
                continue;
            float distance = sqrt((neighbor.position[0] - boid.position[0]) * (neighbor.position[0] - boid.position[0]) + (neighbor.position[1] - boid.position[1]) * (neighbor.position[1] - boid.position[1]) + (neighbor.position[2] - boid.position[2]) * (neighbor.position[2] - boid.position[2]));
            // check that we don't consider ourselves a neighbor
            if (distance > 0) {
                Vector3D difference = neighbor.position - boid.position;
                
                // check whether neighboring boid is within separation distance
                if (distance < separationRadius) {
                    separationVelocity -= (difference / (distance * distance));

                // check whether neighboring boid is within alignment/cohesion distance
                } else if (distance < radius) {
                    averagePosition += neighbor.position;
                    averageVelocity += neighbor.velocity;
                    neighbors++;
                }
            }
        }


        for (Boid &pred : boids) {
            if(boid.isPredator || !pred.isPredator)
                continue;
            Vector3D diff = boid.position - pred.position;
            Vector3D pred_dxyz = Vector3D();
            int pred_count = 0;
            float distance = diff.norm();
            if (distance < predRange) {
                pred_dxyz += diff;
                pred_count++;
            }
            if (pred_count>0) {
                if (pred_dxyz[2] > 0) {
                    boid.acceleration += Vector3D(0, 0, predTurnFactor);
                }
                if (pred_dxyz[2] < 0) {
                    boid.acceleration -= Vector3D(0, 0, predTurnFactor);
                }
                if (pred_dxyz[1] > 0) {
                    boid.acceleration += Vector3D(0, predTurnFactor, 0);
                }
                if (pred_dxyz[1] < 0) {
                    boid.acceleration -= Vector3D(0, predTurnFactor, 0);
                }
                if (pred_dxyz[0] > 0) {
                    boid.acceleration += Vector3D(predTurnFactor, 0, 0);
                }
                if (pred_dxyz[0] < 0) {
                    boid.acceleration -= Vector3D(predTurnFactor, 0, 0);
                }
            }
        }

        if (boid.isPredator) {
            int close_ix = -1;
            double closest_dist = __DBL_MAX__; 

            for (int p = 0; p < boids.size(); p++) {
                if(boids[p].isPredator)
                    continue;
                double dist = (boid.position-boids[p].position).norm();
                if(dist < closest_dist) {
                    closest_dist = dist;
                    close_ix = p;
                }
            }
            boid.acceleration += hungies * (boids[close_ix].position - boid.position);
        }    
        
       
        if(!boid.isPredator){
        // account for cohesion, alignment, and separation
            if (neighbors > 0) {
                averagePosition = averagePosition/neighbors;
                averageVelocity = averageVelocity/neighbors;
                
                boid.acceleration += (cohesionFactor * averagePosition - boid.position) + (alignmentFactor * averageVelocity - boid.velocity);
            }
            boid.acceleration += separationFactor * separationVelocity;
        }
    
        // small random offset to account for stochastic effects
        Vector3D randomForce = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        boid.acceleration += randomForce/100;
        
        // add bias for left or right sides of screen
        for (Boid &boid : right) {
            boid.acceleration += Vector3D((1 - bias) * boid.acceleration[0] + bias, 0, 0);
        }
        for (Boid &boid : left) {
            boid.acceleration += Vector3D((1 - bias) * boid.acceleration[0] - bias, 0, 0);
        }
        
//         check boundaries
        if (boid.position[2] > 0.1) {
            boid.acceleration -= Vector3D(0, 0, boundaryFactor);
        }
        if (boid.position[2] < -0.1) {
            boid.acceleration += Vector3D(0, 0, boundaryFactor);
        }
        if (boid.position[1] > 0.1) {
            boid.acceleration -= Vector3D(0, boundaryFactor, 0);
        }
        if (boid.position[1] < -0.1) {
            boid.acceleration += Vector3D(0, boundaryFactor, 0);
        }
        if (boid.position[0] > 0.1) {
            boid.acceleration -= Vector3D(boundaryFactor, 0, 0);
        }
        if (boid.position[0] < 0.1) {
            boid.acceleration += Vector3D(boundaryFactor, 0, 0);
        }
        
        // update boid position and velocity
        boid.velocity += boid.acceleration * delta_t;
        boid.position += boid.velocity * delta_t;
    
        // adjust speed if out of bounds
        double speed = sqrt(boid.velocity[0] * boid.velocity[0] + boid.velocity[1] * boid.velocity[1] + boid.velocity[2] * boid.velocity[2]);

        if (speed > maxSpeed) {
            boid.velocity.normalize();
            boid.velocity *= maxSpeed;
        }
        if (speed < minSpeed) {
            boid.velocity.normalize();
            boid.velocity *= minSpeed;
        }
        
        // reset boid acceleration
        boid.acceleration = Vector3D(0, 0, 0);
        
        // find new density dependent radius
        double s = 0.5 * delta_t;
        double d = radius - (0.05 * neighbors);
        radius = max(separationRadius, (((1 - s) * radius) + (s * d)));
    }
}

void Flock::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4): Build a spatial map out of all of the point masses.

}

void Flock::self_collide(Boid &boid, double simulation_steps) {
  // TODO (Part 4): Handle self-collision for a given point mass.

}

float Flock::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.

  return 0.f; 
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Flock::reset() {
  Boid *boid = &boids[0];
  for (int i = 0; i < boids.size(); i++) {
    boid->position = boid->start_position;
    boid++;
  }
}
