#Affine Tranformations for `Boost.Astronomy`
#### GSoC 2020 - programming competency test

#### Introduction:

Affine transformations are used to perform tranformations like translate, rotate, shear, scale, etc over vector spaces.
Class `affine_tranformation`  has  been implemented to be accomodated in namespace `boost::astronomy::coordinate`. It can apply tranformations to `spherical`, ` spherical_equatorial`, and `cartesian` representations.

#### Transformations supported :

* scale
* shear
* translate
* rotate_x
* rotate_y
* rotate_z
*Note : Order in which transformations are applied does not matter except for order of rotation about each axis*

After the above tranformations have been used according to neccessity, they need to be applied to the required coordinate using the `transform` function.

#### Usage example:
```c++
	//An object of affine_transformation to store the transformation matrix
	affine_transformation affine; 
	
	//required transformations are applied to the transformation matrix
    affine.rotate_y(30);
    affine.translate(5.6,4.4,3.7);
    affine.scale(4.5,7.6,11.4);
    
    //Creation of points to be transformed
    auto point1 = make_cartesian_representation(3.0 * meters, 5.0 * si::kilo *meters, 4.0 * si::mega * meters);
    auto point2 = make_spherical_representation(45.0 * bud::degrees, 45.0 * bud::degrees, 1.0 * meters);
    auto point3 = make_spherical_equatorial_representation(3.0 * bud::degrees, 50.0 * bud::degrees, 40.0 * meters);

	//tranformed points can be obtained through transform function
    auto point4 = affine.transform(point1);
    auto point5 = affine.transform(point2);
    auto point6 = affine.transform(point3);
``` 
