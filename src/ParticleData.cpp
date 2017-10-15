#include "ParticleData.hpp"




// ============================================================================
ParticleData::ParticleData()
{
    int N = 10000;

    for (int n = 0; n < N; ++n)
    {
        auto p = Particle();
        p.position = double(n) / (N - 1);
        p.velocity = 0;
        p.momentum = 0;
        particles.push_back(p);
    }
}
