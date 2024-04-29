#ifndef FLOCK_H
#define FLOCK_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "boid.h"
#include "CGL/CGL.h"
#include "CGL/misc.h"
//#include "flockMesh.h"
#include "collision/collisionObject.h"

using namespace CGL;
using namespace std;

struct FlockParameters {
    FlockParameters() {}
    FlockParameters(double radius, double separationRadius, double boundaryFactor, double separationFactor, double alignmentFactor, double cohesionFactor, double maxSpeed, double minSpeed, bool enable_separation, bool enable_alignment, bool enable_cohesion)
      : radius(radius), separationRadius(separationRadius), boundaryFactor(boundaryFactor), separationFactor(separationFactor), alignmentFactor(alignmentFactor), cohesionFactor(cohesionFactor), maxSpeed(maxSpeed), minSpeed(minSpeed), enable_separation(enable_separation), enable_alignment(enable_alignment), enable_cohesion(enable_cohesion) {}
  ~FlockParameters() {}
    
    // Global simulation parameters

    bool enable_separation = true;
    bool enable_alignment = true;
    bool enable_cohesion = true;

    double radius = 2;
    double separationRadius = 0.1;
    double boundaryFactor = 0.4;
    double separationFactor = 0.05;
    double alignmentFactor = 0.005;
    double cohesionFactor = 0.05;
    double predRange = 50;
    double predTurnFactor = 0.15;
    double maxSpeed = 1.6;
    double minSpeed = 0.8;
    double hungies = 1.0;
    double windPower = 0;
    int num_boids = 200;
};

struct Flock {
    Flock() {}
    Flock(int num_boids);
  ~Flock();

  void buildGrid();
  bool fieldOfView(Vector3D behind, Vector3D difference, double angle);
  void simulate(double frames_per_sec, double simulation_steps, FlockParameters *cp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects);

  void reset();
  void buildFlockMesh();

  void build_spatial_map();
  void self_collide(Boid &boid, double simulation_steps);
  float hash_position(Vector3D pos);
  void update(Boid &boid);

  // Flock properties
    int num_boids = 200;
    int sim_step = 0;

  // Flock components
    vector<Boid> boids;
//    FlockMesh *flockMesh;

  // Spatial hashing
  unordered_map<float, vector<Boid *> *> map;
};

#endif /* FLOCK_H */
