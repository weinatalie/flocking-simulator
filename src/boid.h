#ifndef BOID_H
#define BOID_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"

using namespace CGL;

struct Boid {
  Boid(Vector3D position, Vector3D velocity, Vector3D acceleration)
      : start_position(position), position(position), velocity(velocity), acceleration(acceleration)
    {}

  Vector3D normal();

  // static values
  Vector3D start_position;

  // dynamic values
  Vector3D position;
  Vector3D velocity;
  Vector3D acceleration;
};

#endif /* BOID_H */
