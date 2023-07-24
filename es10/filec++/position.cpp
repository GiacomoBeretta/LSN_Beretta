#include "position.hpp"

Position ::Position() {
  x_ = 0;
  y_ = 0;
}

Position ::Position(const double x, const double y) {
  x_ = x;
  y_ = y;
}

Position ::Position(const Position &p) {
  x_ = p.getX();
  y_ = p.getY();
}

void Position ::setX(const double x) { x_ = x; }

void Position ::setY(const double y) { y_ = y; }

double Position ::distance(const Position &p) const {
  double dx = p.getX() - x_;
  double dy = p.getY() - y_;
  return sqrt(dx * dx + dy * dy);
}

double Position::distanceOnSphere(const Position &p) const {
  double phi1, theta1, x1, y1, z1, phi2, theta2, x2, y2, z2, dx, dy, dz, d;
  phi1 = getX();
  theta1 = getY();
  x1 = (sin(theta1)) * cos(phi1);
  y1 = (sin(theta1)) * sin(phi1);
  z1 = cos(theta1);

  phi2 = p.getX();
  theta2 = p.getY();
  x2 = (sin(theta2)) * cos(phi2);
  y2 = (sin(theta2)) * sin(phi2);
  z2 = cos(theta2);

  dx = x1 - x2;
  dy = y1 - y2;
  dz = z1 - z2;
  d = sqrt(dx * dx + dy * dy + dz * dz);
  return 2.0 * asin(d / 2.0);
}

double Position ::getRx(const Position &p) const { return p.getX() - x_; }

double Position ::getRy(const Position &p) const { return p.getY() - y_; }

bool Position ::operator==(const Position &p) const {
  if (x_ == p.getX() && y_ == p.getY())
    return true;
  else
    return false;
}

bool Position ::operator!=(const Position &p) const { return not(*this == p); }

Position &Position ::operator=(const Position &p) {
  if (this != &p) {
    x_ = p.getX();
    y_ = p.getY();
  }

  return *this;
}
