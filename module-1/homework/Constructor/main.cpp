#include <stdio.h>
#include <stdint.h>

class Bullet
{
public:
    Bullet() = delete;

    Bullet(int32_t x, int32_t y, int32_t z)
      : m_x(x), m_y(y), m_z(z)
    { };

    Bullet(const Bullet& bullet)
      : m_x(bullet.m_x), m_y(bullet.m_y), m_z(bullet.m_z)
    { 

    };

    ~Bullet() = default;

private:
    int32_t m_x;

    int32_t m_y;

    int32_t m_z;
};

int main()
{
    int32_t x = 0;
    int32_t y = 1;
    int32_t z = 3;

    Bullet bullet(x,y,z);

    Bullet bullet1 = bullet;

    return 0;
}