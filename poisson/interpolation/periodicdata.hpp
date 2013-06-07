// Jan Becker
// 2012-07-10

#pragma once

/**
 * Data needed for mapping a periodic velocity / height profile to a time
 * dependent one.
 *
 * Used for Periodic Balakotaiah + ALE and DNS + ALE
 */
struct PeriodicData {

  double c; //phase-velocity
  double deltaX;//increment x-direction
  double deltaY;//increment y-direction
  int NY;//#gridpoints y-direction
  int NX;//#gridpoints x-direction

  PeriodicData()
      : c(0.385), deltaX(0.000209), deltaY(0.0001), NY(41), NX(102)
  {}

};
