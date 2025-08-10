/*
	Classes which define rotations about the x, y, and z axes. Using these,
	any arbitrary orientation can be described. Methods implemented for vector multiplication
	as well as generating the inverse of the rotation.
*/
#include "Rotation.h"

namespace Mask {

XRotation::XRotation() :
m_angle(0)
{
	GenerateMatrix();
}

XRotation::XRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

XRotation::~XRotation() {}

void XRotation::GenerateMatrix() {
	m_matrix[0][0] = 1.0; m_matrix[0][1] = 0.0; m_matrix[0][2] = 0.0;
	m_matrix[1][0] = 0.0; m_matrix[1][1] = std::cos(m_angle); m_matrix[1][2] = -std::sin(m_angle);
	m_matrix[2][0] = 0.0; m_matrix[2][1] = std::sin(m_angle); m_matrix[2][2] = std::cos(m_angle);
}

YRotation::YRotation() :
m_angle(0)
{
	GenerateMatrix();
}

YRotation::YRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

YRotation::~YRotation() {}

void YRotation::GenerateMatrix() {
	m_matrix[0][0] = std::cos(m_angle); m_matrix[0][1] = 0.0; m_matrix[0][2] = -std::sin(m_angle);
	m_matrix[1][0] = 0.0; m_matrix[1][1] = 1.0; m_matrix[1][2] = 0.0;
	m_matrix[2][0] = std::sin(m_angle); m_matrix[2][1] = 0.0; m_matrix[2][2] = std::cos(m_angle);
}


ZRotation::ZRotation() :
m_angle(0)
{
	GenerateMatrix();
}

ZRotation::ZRotation(double angle) :
m_angle(angle)
{
	GenerateMatrix();
}

ZRotation::~ZRotation() {}

void ZRotation::GenerateMatrix() {
	m_matrix[0][0] = std::cos(m_angle); m_matrix[0][1] = -std::sin(m_angle); m_matrix[0][2] = 0.0;
	m_matrix[1][0] = std::sin(m_angle); m_matrix[1][1] = std::cos(m_angle); m_matrix[1][2] = 0.0;
	m_matrix[2][0] = 0.0; m_matrix[2][1] = 0.0; m_matrix[2][2] = 1.0;
}

};