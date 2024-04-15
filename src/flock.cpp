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
    for (int i = 0; i < num_boids; i++)
    {
        Vector3D position = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        Vector3D velocity = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        Vector3D acceleration = Vector3D(0, 0, 0);
        Boid boid = Boid(position, velocity, acceleration);
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
    double maxSpeed;
    double minSpeed;
    double bias;
    bool cohesion;
    bool alignment;
    bool separation;
    
    radius = fp->radius;
    separationRadius = fp->separationRadius;
    boundaryFactor = fp->boundaryFactor;
    separationFactor = fp->separationFactor;
    alignmentFactor = fp->alignmentFactor;
    cohesionFactor = fp->cohesionFactor;
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
    
    // iterate through all boids
    for (Boid &boid : boids) {
        int neighbors = 0;
        Vector3D averagePosition = Vector3D(0, 0, 0);
        Vector3D averageVelocity = Vector3D(0, 0, 0);
        Vector3D separationVelocity = Vector3D(0, 0, 0);
        
        // iterate through all neighboring boids
        for (Boid &neighbor : boids) {
            
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
        
        // account for cohesion, alignment, and separation
        if (neighbors > 0) {
            averagePosition = averagePosition/neighbors;
            averageVelocity = averageVelocity/neighbors;
            
            boid.acceleration += (cohesionFactor * averagePosition - boid.position) + (alignmentFactor * averageVelocity - boid.velocity);
        }
        boid.acceleration += separationFactor * separationVelocity;
        
        // random offset to account for stochastic effects
        Vector3D randomForce = Vector3D((rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1, (rand() / double(RAND_MAX)) * 2 - 1);
        boid.acceleration += randomForce/100;
        
        // update boid position and velocity
        boid.position += boid.velocity * delta_t;
        boid.velocity += boid.acceleration * delta_t;
        
        // check boundaries
        if (boid.position[2] >= 1) {
            boid.velocity -= Vector3D(0, 0, boundaryFactor);
        }
        if (boid.position[2] <= -1) {
            boid.velocity += Vector3D(0, 0, boundaryFactor);
        }
        if (boid.position[1] >= 1) {
            boid.velocity -= Vector3D(0, boundaryFactor, 0);
        }
        if (boid.position[1] <= -1) {
            boid.velocity += Vector3D(0, boundaryFactor, 0);
        }
        if (boid.position[0] >= 1) {
            boid.velocity -= Vector3D(boundaryFactor, 0, 0);
        }
        if (boid.position[0] <= -1) {
            boid.velocity += Vector3D(boundaryFactor, 0, 0);
        }
    
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

//void Flock::buildFlockMesh() {
//  if (boids.size() == 0) return;
//
//  FlockMesh *flockMesh = new FlockMesh();
//  vector<Triangle *> triangles;
//
//  // Create vector of triangles
//  for (int y = 0; y < num_height_points - 1; y++) {
//    for (int x = 0; x < num_width_points - 1; x++) {
//      PointMass *pm = &point_masses[y * num_width_points + x];
//      // Get neighboring point masses:
//      /*                      *
//       * pm_A -------- pm_B   *
//       *             /        *
//       *  |         /   |     *
//       *  |        /    |     *
//       *  |       /     |     *
//       *  |      /      |     *
//       *  |     /       |     *
//       *  |    /        |     *
//       *      /               *
//       * pm_C -------- pm_D   *
//       *                      *
//       */
//
//      float u_min = x;
//      u_min /= num_width_points - 1;
//      float u_max = x + 1;
//      u_max /= num_width_points - 1;
//      float v_min = y;
//      v_min /= num_height_points - 1;
//      float v_max = y + 1;
//      v_max /= num_height_points - 1;
//
//      PointMass *pm_A = pm                       ;
//      PointMass *pm_B = pm                    + 1;
//      PointMass *pm_C = pm + num_width_points    ;
//      PointMass *pm_D = pm + num_width_points + 1;
//
//      Vector3D uv_A = Vector3D(u_min, v_min, 0);
//      Vector3D uv_B = Vector3D(u_max, v_min, 0);
//      Vector3D uv_C = Vector3D(u_min, v_max, 0);
//      Vector3D uv_D = Vector3D(u_max, v_max, 0);
//
//
//      // Both triangles defined by vertices in counter-clockwise orientation
//      triangles.push_back(new Triangle(pm_A, pm_C, pm_B,
//                                       uv_A, uv_C, uv_B));
//      triangles.push_back(new Triangle(pm_B, pm_C, pm_D,
//                                       uv_B, uv_C, uv_D));
//    }
//  }
//
//  // For each triangle in row-order, create 3 edges and 3 internal halfedges
//  for (int i = 0; i < triangles.size(); i++) {
//    Triangle *t = triangles[i];
//
//    // Allocate new halfedges on heap
//    Halfedge *h1 = new Halfedge();
//    Halfedge *h2 = new Halfedge();
//    Halfedge *h3 = new Halfedge();
//
//    // Allocate new edges on heap
//    Edge *e1 = new Edge();
//    Edge *e2 = new Edge();
//    Edge *e3 = new Edge();
//
//    // Assign a halfedge pointer to the triangle
//    t->halfedge = h1;
//
//    // Assign halfedge pointers to point masses
//    t->pm1->halfedge = h1;
//    t->pm2->halfedge = h2;
//    t->pm3->halfedge = h3;
//
//    // Update all halfedge pointers
//    h1->edge = e1;
//    h1->next = h2;
//    h1->pm = t->pm1;
//    h1->triangle = t;
//
//    h2->edge = e2;
//    h2->next = h3;
//    h2->pm = t->pm2;
//    h2->triangle = t;
//
//    h3->edge = e3;
//    h3->next = h1;
//    h3->pm = t->pm3;
//    h3->triangle = t;
//  }
//
//  // Go back through the cloth mesh and link triangles together using halfedge
//  // twin pointers
//
//  // Convenient variables for math
//  int num_height_tris = (num_height_points - 1) * 2;
//  int num_width_tris = (num_width_points - 1) * 2;
//
//  bool topLeft = true;
//  for (int i = 0; i < triangles.size(); i++) {
//    Triangle *t = triangles[i];
//
//    if (topLeft) {
//      // Get left triangle, if it exists
//      if (i % num_width_tris != 0) { // Not a left-most triangle
//        Triangle *temp = triangles[i - 1];
//        t->pm1->halfedge->twin = temp->pm3->halfedge;
//      } else {
//        t->pm1->halfedge->twin = nullptr;
//      }
//
//      // Get triangle above, if it exists
//      if (i >= num_width_tris) { // Not a top-most triangle
//        Triangle *temp = triangles[i - num_width_tris + 1];
//        t->pm3->halfedge->twin = temp->pm2->halfedge;
//      } else {
//        t->pm3->halfedge->twin = nullptr;
//      }
//
//      // Get triangle to bottom right; guaranteed to exist
//      Triangle *temp = triangles[i + 1];
//      t->pm2->halfedge->twin = temp->pm1->halfedge;
//    } else {
//      // Get right triangle, if it exists
//      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
//        Triangle *temp = triangles[i + 1];
//        t->pm3->halfedge->twin = temp->pm1->halfedge;
//      } else {
//        t->pm3->halfedge->twin = nullptr;
//      }
//
//      // Get triangle below, if it exists
//      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
//        Triangle *temp = triangles[i + num_width_tris - 1];
//        t->pm2->halfedge->twin = temp->pm3->halfedge;
//      } else {
//        t->pm2->halfedge->twin = nullptr;
//      }
//
//      // Get triangle to top left; guaranteed to exist
//      Triangle *temp = triangles[i - 1];
//      t->pm1->halfedge->twin = temp->pm2->halfedge;
//    }
//
//    topLeft = !topLeft;
//  }
//
//  flockMesh->triangles = triangles;
//  this->flockMesh = flockMesh;
//}
