#ifndef __position_hpp__
#define __position_hpp__

#include <cmath>
#include <iostream>

const double pi = 3.14159265358979323846;

class Position {

public:
  // costruttori
  Position();
  Position(const double, const double);
  Position(const Position &);
  // non serve un distruttore basta quello del compilatore

  double getX() const { return x_; } // Coordinate cartesiane
  double getY() const { return y_; }
  void setX(const double);
  void setY(const double);
  double distance(const Position &) const; // distanza da un altro punto
  double distanceOnSphere(
      const Position &) const; // distanza di due punti su una sfera di
                               // raggio R calcolata come geodetica

  // vettore differenza r2-r1
  double getRx(const Position &) const;
  double getRy(const Position &) const;
  // double getRz(const Position &) const;

  bool operator==(const Position &) const;
  bool operator!=(const Position &) const;
  Position &operator=(const Position &);

protected:
  double x_, y_;
};

#endif // __position_hpp__
