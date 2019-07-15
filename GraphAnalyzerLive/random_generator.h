#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

class random_generator
{
private:
    long long rnd;
    long long range;
    long long low_range;
    long long high_range;
    long long low;
    long long high;
    long long remainder;
    double temp_gauss;
    double save_gauss;
public:
    random_generator();
    long long random(long l, long h);
    double random_double();
    double random_gauss();
};

#endif // RANDOM_GENERATOR_H
