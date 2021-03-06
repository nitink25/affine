#ifndef BOOST_ASTRONOMY_COORDINATE_AFFINE_TRANSFORMATION_HPP
#define BOOST_ASTRONOMY_COORDINATE_AFFINE_TRANSFORMATION_HPP

//!Header to facilitate basic affine transforms like shear, scale, rotate, translate.
//Supports Cartesian, Spherical and Spherical Equatorial Represenations.

#include <iostream>
#include <cmath>
#include <array>

#include <boost/astronomy/coordinate/cartesian_representation.hpp>
#include <boost/astronomy/coordinate/spherical_representation.hpp>
#include <boost/astronomy/coordinate/spherical_equatorial_representation.hpp>

namespace boost { namespace astronomy { namespace coordinate {

class affine_transformation
{
private:    
   
    //!tranformation matrix
    std::array<std::array<long double,4>, 4> transMatrix;

public:

    //!Constructor to create an object and to initialize transMatrix
    affine_transformation()
    {
        //tranformation matrix initially declared as identity
        this->transMatrix = {{ {1,0,0,0},
                               {0,1,0,0},
                               {0,0,1,0},
                               {0,0,0,1} }};
    }

    //!function to update the transMatrix with a new transform
    void update( std::array<std::array<long double,4>, 4> const& matrix )
    {
        std::array<std::array<long double,4>, 4> temp;

        //matrix multiplication
        for( int i = 0; i < 4; i++ )
        {
            for( int j = 0; j < 4; j++ )
            {
                temp.at(i).at(j) = 0;
                for ( int k = 0; k < 4; k++ )
                {
                    temp.at(i).at(j) += this->transMatrix.at(i).at(k) * matrix.at(k).at(j);
                }
            }
        }

        this->transMatrix = temp;
    }

    //!function to apply scaling
    void scale( long double sx, long double sy, long double sz )
    {
        std::array<std::array<long double,4>, 4> temp = {{ {sx, 0 , 0 , 0},
                                                           {0 , sy, 0 , 0},
                                                           {0 , 0 , sz, 0},
                                                           {0 , 0 , 0 , 1} }};

        this->update(temp);
    }

    //!function to apply shear
    void shear( long double hxy, long double hxz,
                long double hyx, long double hyz, 
                long double hzx, long double hzy )
    {
        std::array<std::array<long double,4>, 4> temp = {{ { 1 ,hxy,hxz,0},
                                                           {hyx, 1 ,hyz,0},
                                                           {hzx,hzy, 1 ,0},
                                                           { 0 , 0 , 0 ,1} }};

        this->update(temp);
    }

    //!function to apply rotation about x-axis
    //only degrees accepted for all rotation functions
    void rotate_x( long double x )
    {
        std::array<std::array<long double,4>, 4> temp = {{ {1 , 0 , 0 , 0},
                                                           {0 ,std::cos(x),-std::sin(x) ,0},
                                                           {0 ,std::sin(x), std::cos(x),0},
                                                           {0 , 0 , 0 , 1} }};

        this->update(temp);
    }

    //!function to apply rotation about y-axis
    void rotate_y( long double y )
    {
        std::array<std::array<long double,4>, 4> temp = {{ { std::cos(y), 0 ,std::sin(y), 0},
                                                           {0 , 1 , 0 , 0},
                                                           {-std::sin(y), 0 ,std::cos(y),0},
                                                           {0 , 0 , 0 , 1} }};

        this->update(temp);
    }

    //!function to apply rotation about z-axis
    void rotate_z( long double z )
    {
        std::array<std::array<long double,4>, 4> temp = {{ {std::cos(z),-std::sin(z), 0 , 0},
                                                           {std::sin(z), std::cos(z), 0 , 0},
                                                           {0 , 0 , 1 , 0},
                                                           {0 , 0 , 0 , 1} }};

        this->update(temp);
    }

    //!function to apply translation
    void translate( long double dx, long double dy, long double dz )
    {
        std::array<std::array<long double,4>, 4> temp = {{ {1 ,0 ,0 ,dx},
                                                           {0 ,1 ,0 ,dy},
                                                           {0 ,0 ,1 ,dz},
                                                           {0 ,0 ,0 ,1 } }};

        this->update(temp);
    }

    //!function to print the transformation matrix
    void display()
    {
        for( int i = 0 ; i < 4 ; i++ )
        {
            for( int j = 0 ; j < 4 ; j++ )
            {
                std::cout << this->transMatrix.at(i).at(j) << " ";
            }

            std::cout<<"\n";   
        } 
    }

    //!function to apply the transformations
    template
    <
        template<typename ...> class Representation,
        typename ...Args
    >
    auto transform( Representation<Args...> const& vector )
    {
        auto tempVector = make_cartesian_representation(vector);
        typedef decltype(tempVector) temp_type;

        //!homogeneous coordinates matrix
        std::array<long double, 4> homoMatrix = { tempVector.get_x().value(),
                                                  tempVector.get_y().value(),
                                                  tempVector.get_z().value(),
                                                  1 };

        std::array<long double, 4> tempMatrix = {0,0,0,0};

        for ( int i = 0 ; i < 4 ; i++ )
        {
            for ( int j = 0 ; j < 4 ; j++ )
            {    
                tempMatrix.at(i) += this->transMatrix.at(i).at(j) * homoMatrix.at(j);
            }   
        }

        tempVector.set_x( static_cast<typename temp_type::type>(tempMatrix.at(0)) * 
            typename temp_type::quantity1::unit_type() );
        tempVector.set_y( static_cast<typename temp_type::type>(tempMatrix.at(1)) * 
            typename temp_type::quantity2::unit_type() );
        tempVector.set_z( static_cast<typename temp_type::type>(tempMatrix.at(2)) * 
            typename temp_type::quantity3::unit_type() );

        return Representation<Args...>(tempVector);
    }

};

}}} //namespace boost::astronomy::coordinate
#endif // !BOOST_ASTRONOMY_COORDINATE_AFFINE_TRANSFORMATION_HPP