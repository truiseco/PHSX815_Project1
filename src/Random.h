/*************************************************************
* @author   Triston Ruiseco
* @file     Random.h
* @date     02/08/2021
* @brief    Class to generate pseudo-random numers.
*************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <cmath>

class Random
{
  public:

    /**
    * @pre none
    * @post random object instantiated
    * @param seed: optional starting point for number randomization
    */
    Random(unsigned long long seed = 5555);

    /**
    * @pre none
    * @post none
    * @return random 64 bit int
    */
    unsigned long long int64();

    /**
    * @pre none
    * @post none
    * @return random 32 bit int
    */
    unsigned int32();

    /**
    * @pre none
    * @post none
    * @return uniform random double between 0 and 1
    */
    double rand();

    /**
    * @pre none
    * @post none
    * @return random integer 0,1 according to Bernoulli dist with prob p
    * @param p: probability of returning 1
    */
    int Bernoulli(double p = 0.5);

    /**
    * @pre none
    * @post none
    * @return random double 0-infinity according to an exponential distribution
    * @param beta: expected value of distribution, also characterizing parameter
                   for sharpening/softening of distribution
    */
    double Exponential(double beta = 1.0);

    /**
    * @pre none
    * @post none
    * @return random integer 0,1,..,n according to uniform categorical
                     distribution with n categories
    * @param n: categories in distribution
    */
    int Categorical(int n = 6);

    /**
    * @pre none
    * @post none
    */
    virtual ~Random(){};

  private:

    /*
    * @brief private members to give memory to the RNG sequence
    */
    unsigned long long m_u;
    unsigned long long m_v;
    unsigned long long m_w;
};
#endif
