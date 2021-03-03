/*************************************************************
* @author   Triston Ruiseco
* @file     Random.cpp
* @date     02/08/2021
* @brief    Implementation file for Random.h.
*************************************************************/
#include "Random.h"

Random::Random(unsigned long long seed /*=5555*/){
  m_v = 4101842887655102017LL;
  m_w = 1;

  m_u = seed ^ m_v;
  int64();
  m_v = m_u;
  int64();
  m_w = m_v;
  int64();
}

unsigned long long Random::int64(){
  m_u = m_u * 2862933555777941757LL + 7046029254386353087LL;
  m_v ^= m_v >> 17;
  m_v ^= m_v << 31;
  m_v ^= m_v >> 8;
  m_w = 4294957665U*(m_w & 0xffffffff) + (m_w >> 32);
  unsigned long long x = m_u ^ (m_u << 21);
  x ^= x >> 35;
  x ^= x << 4;
  return (x + m_v) ^ m_w;
}

unsigned Random::int32(){
  return (unsigned int)int64();
}

double Random::rand(){
  return 5.42101086242752217E-20 * int64();
}

int Random::Bernoulli(double p){
  if(p < 0. || p > 1.)
    return 1;

  if(rand() < p)
    return 1;
  else
    return 0;
}

double Random::Exponential(double rate){
  if(rate <= 0.)
    rate = 1.;

  double R = rand();
  while(R <= 0.)
    R = rand();

  double X = -log(R)/rate;

  return X;
}

int Random::Categorical(int n){
  if(n < 3)
    n = 3;

  double R = 1.0 + n*rand();
  return (int)R;
}
