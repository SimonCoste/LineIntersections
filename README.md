# LineIntersections



This package contains resources for the generation of 2-dimensional random line processes, as well as their intersection points. In particular, there is a tailored version of the Bentley-Ottman algorithm.  

![Poisson intersection process](docs/poisson_intersection.png)

## Line objects

A two-dimensional line is uniquely parametrized by the coordinates of the projection of the origin on the line, that is, by two numbers (theta, r) where theta is an angle and r a nonnegative number. This is the content of the `Line` type.

For the moment there are two random line generators : Poisson, and triangular. 

*Example*

```julia
radius = 1 ; number_of_lines = 50
H = hyperplanes_poisson(radius, number_of_lines)
intersection_points = intersections_bentley_ottman(H)
```

## Intersections

The function intersections_bentley_ottman takes as input an array of lines and an optional radius R = 1 ; the output is a $(2, m)$ array where each column is the coordinate of one of the m intersection points between the lines. The intersection points are sorted by increasing x-coordinate. 

## Todo : 
- Bentley-Ottman algo is flawed, it finds all intersection points but overcounts some of them. 
- Vertical lines
- Multiple intersections
- Implement other line processes
- Cool animated pic


