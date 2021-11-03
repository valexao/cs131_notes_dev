---
title: Camera Parameters
keywords: camera parameters, triangulation, multi-view geometry, epipolar geometry
order: 11 # Lecture number for 2021
---

# CS 131 Lecture Notes
### October 26, 2021
#### Valexa Orelien, Alexis Lowber, Jacob Castaneda, Colin Norick, Adonis Pugh, Jason Chen

# Multi-view Geometry

In these notes, we will answer the question: what are the geometric constraints governing multiple views of the same 3D world?

## Why multi-view geometry?

When there are *known* cameras that output images with estimated 2D correspondences, we can find 3D points and triangulate them to create a ***recovering structure***. 

We can see this through a model below, where we know that the camera intrisic and extrinsic $R$ and $t$ are the rotation and translation of each camera relative to each other.

![](https://i.imgur.com/hRCSAkc.png | width=100)


Recall the definition of ***motion*** as estimating the $R$ and $t$ parameters for a set of cameras given 2D correspondences to find where they are in the world.

There are certain geometric constraints that allow us to calculate the mathematical relationship between recovering structure and motion. Triangulation, epipolar geometry, and stereo are all concepts that help to describe this relationship. In this document, we will be going into further detail about triangulation and epipolar geometry.

## Triangulation
### Linear Approach
We have a problem where we are given projections of 3D points in two or more images with known matrices and need to use that information to find the coordinates of a point. There is no way to be absolutely sure about everything about the images  because we don't necessarily know about absolute position, orientation, or scale. This ambiguity is observed because it is possible to apply a projective transformation that is represented by a 4 × 4 matrix $H$ to each point $x_i$, and on the right of each camera matrix $P_j$ , without changing the projected image points such that $P_jx_i = (P_j H^{−1} )(H*x_i)$ (Hartley and Zisserman).

Let’s say we have two images, where we know the 3D points $x$ and $x'$. Let $X$ represent the coordinates of the point that we are trying to find, and we can use this model to visually represent how we are trying to find $X$.

![Model with two cameras](https://i.imgur.com/KPlfgTv.png)

Recall that we are using the pinhole camera model, so we are using a camera sensor, which is a ray that represents going from the camera center ($O$ for the first camera and $O'$ for the second camera in this case) to the pixel. This ray is like a linear projection that allows us to project from 2 dimensional to 3 dimensional points using a pinhole camera projection matrix, also known as a $P$ matrix. This will tell us that there is the 3D point $X$ somewhere along that line at an unknown depth. We would apply this to both cameras in this case - if there are more 2 cameras then the camera sensor should be calculated for each camera - and trace that ray for the camera, which again tells us the direction of the 3D point $X$ at an unknown depth. 

The next step might be to ask, where do these rays intersect? If we know that they correspond to the same three points, then we can infer that they should intersect at a 3D point. That actually is the definition of ***triangulation*** - that when are given two or more points that correspond to each other, the rays should intersect at a 3D point. 

However, the likelihood of these rays intersecting can be very slim if we are not in a very idealized case. In practice, this can be very difficult due to noise and distance - the further cameras are away from a certain point, the *more sensitive they are to small variations* which makes it harder to find a ray that satisfies our criteria and intersects with $X$, especially given that we are estimating $x$.

One way to resolve this and “create” idealized conditions is by approximating idealized points and in turn approximating “perfect” rays. We can find the shortest segment where those two rays intersect and let $X$ be the midpoint of that segment in order to generate an estimate. This is a *geometric approach*, where we take the average of the points and estimate that they represent the ideal 3D point.

As aforementioned, there is a major flaw with the geometric approach to triangulation:
Uncertainty, or noise, varies with distance. That is to say, an estimation for a point that is farther from the image planes is more likely to be inaccurate. In statistics, this condition is known as ***heteroscedasticity***.

We can formalize the geometric approach using linear algebra:

First, we know that the optical center ($O$), the 2D point on the image plane ($x$), and the real 3D point($X$) lie on the same ray. They are therefore colinear. So, we can relate these values using a projection matrix $P$ and a constant $\lambda$ to account for scale):


\begin{equation}
    \lambda x = PX \\
    \lambda ' x' = P'X\\
\end{equation}


Since $x$ and $PX$ (and $x'$ and $P'X$ are colinear, so their cross product is 0):

\begin{equation}
    x \times PX = 0 \\
    x' \times P'X = 0 \\
    [x]_{\times} PX = 0 \\
    [x']_{\times} P'X = 0 
\end{equation}

(See [1] for a quick refresher on cross products)

Since $x = x_1 + x_2 + x_3$ is a 3x1 vector, we have three equations. However, simplifiication of the system gives us only two independent equations for $x$ (and two for $x'$). We can write this in matrix form:

![](https://i.imgur.com/dRhcyh0.png)


We want to find a solution for this matrix, but we need to recognize that the solution will not be perfect due to noise. So, instead of setting the error 0, we will instead use a min function which allows us to get the error (the norm of $AX^2$) as close to zero as possible under the constraint that the norm of $X$ is equal to 1:

![](https://i.imgur.com/m1iVl0N.png)


This is a total least squares problem which can be solved using singular value decomposition.


### Non-Linear Approach

Now that we’ve discussed the linear/algebraic approach to triangulation, we can consider the nonlinear approach. With this approach, we begin by finding our projected estimate of the 3D point using our two rays. Then, our goal becomes finding the $X$ that minimizes the 2D reprojection errors:

![](https://i.imgur.com/hzjnXiE.png)

(Note that mathematical techniques used for this model are very similar to the ones used to solve our linear model.)

### [1] Cross-Product Review
#### Geometric representation of cross product:

![](https://upload.wikimedia.org/wikipedia/commons/4/4e/Cross_product_parallelogram.svg)

As a refresher, the cross product of vectors $a$ and $b$ can be represented by the area of the parallelogram formed by $a$ and $b$. So for collinear vectors, this area is 0.

#### Matrix representation of cross product:
![reference link](https://wikimedia.org/api/rest_v1/media/math/render/png/51136b4d7561e04a62ccec43b72e5162f9341633)

This multiplication results in:

\begin{align}
    (a \times b)_1 = 0*b_1 - a_3*b_2 + a_2*b_3 = - a_3*b_2 + a_2*b_3 \\
    (a \times b)_2 = a_3*b_1 + 0*b_2 - a_1*b_3 = a_3*b_1 - a_1*b_3 \\
    (a \times b)_3 = -a_2*b_1 + a_1*b_2 + 0*b_3 = -a_2*b_1 + a_1*b_2
\end{align}

Now, since $a$ and $b$ are collinear, we know:

\begin{align}
    \frac{a_1}{b_1} = \frac{a_2}{b_2} = \frac{a_3}{b_3}
\end{align}

By cross multiplying the ratios, we get:

\begin{align}
    a_3*b_2 = a_2*b_3 \\
    a_3*b_1 = a_1*b_3 \\
    a_2*b_1 = a_1*b_2 \\
\end{align}

Therefore, 

\begin{align}
a \times b = 
\begin{bmatrix}
0 & 0 & 0
\end{bmatrix}
\end{align}


## Epipolar Geometry

### Epipolar Geometry Setup (Jason)
Epipolar geometry involves the geometry between the three-dimensional scene points, the two cameras, the projection of these points onto the image planes, and underlying observations of the above. The standard epipolar geometry setup is diagrammed below, and the two cameras are positioned at points C and C’.

**Epipoles.** The epipoles are where the projection of both cameras onto the image of the other lie. The projection of the right camera onto the image of the left camera is indicated by the point e, while the projection of the left camera onto the image of the right camera is indicated by the point e’. The epipoles are unique for each stereo system, and they form what is known as a baseline onto which both planes lie. 

**Epipolar Plane.** Now, if you consider the three-dimensional point M, which is known as the scene point, and the two points C and C’ where the cameras lie, you will find that these three points form the epipolar plane. Each scene point 

**Epipolar Lines.** The epipolar lines are the intersections of the epipolar plane, defined above, with the two image planes corresponding to the two cameras. The epipolar lines intersect the plane at the respective epipoles of the two cameras. 

![](https://i.imgur.com/SIl2e7R.png)

#### Example Configurations
*Motion Parallel to an Image Plane*
If we consider the diagram below, we see what happens when the two image planes are parallel to one another. For example, this occurs if you are taking a picture of an object and move directly parallel horizontally and take another picture. This causes the two images to share the same plane. As a result, the epipoles designated by e and e’ are located at the point infinity. This is because the baseline that joins the two camera points is parallel to the axis of the two image planes. Tangentially, the two epipolar lines are also parallel to the axis of the two image planes and are horizontal. 

![](https://i.imgur.com/ICO2BDS.png)

*Example Configuration: Motion Perpendicular to an Image Plane*
Motion perpendicular to an image plane occurs when you take a picture, take a step forward, and take another picture, which is common in cases such as self-driving cars. Thus the baseline is perpendicular to the image plane. By construction, the epipoles are the focus of expansion and move outwards, moving away from the principle point which they coincide with.


### Epipolar Constraint
Let’s consider the possibility that you are concerned with finding the corresponding point on the right image of a point on the left image. If you consider two images below and a single point on the left image xₗ, you will see that the potential matches must lie on a epipolar line formed on the right image. This is because the epipolar plane formed by the points of the two cameras and the already-realized point on the left image intersects with the image on the right at a single line. The search for correspondence is thus concentrated on a single epipolar line rather than the entire right-hand image. 
	Consider two rays projected from the point of the camera onto the two points located on matching epipolar lines. You will see that they may intersect in the same three dimensional space, meaning that they both could be projections of the same three dimensional point X. However, you must still be aware that even if both points satisfy the epipolar constraint of the other, they may not be projections of the same three-dimensional point in the image. 

![](https://i.imgur.com/Zhm4cWA.png)

### Math for epipolar constant: Calibrated Case

For the calibrated case, we can make the assumption that the intrinsic and extrinsic parameters are known for both of the cameras, and by convention, the world is centered at the first camera. Thus the projection matrices for the first camera is equal to P below, where the extrinsic matrix is K, which is intrinsic to the camera itself, and since there is no rotation from the center, the rotation matrix is the identity, and the translation matrix is zero, since there is no translation occuring at the origin. However, for the projection matrices for the second camera, indicated by P', the coordinate relies an extrinsic K' and a rotation matrices R and a translation t from the first camera.
\begin{align}
    P = K[I|0] \\
    P' = K'[R|t]
\end{align}

We can then normalize the projection matrices by pre-multiplying them by the inverse calibration matrices. The ray from O, the camera, to the point x is given by the inverse of the extrinsic matrix of the first camera multiplied by the world coordinate X we are examining, and the same is true for O', where the point x' is given by the inverse of the extrinsic matrix of the other camera multiplied by X.
\begin{align}
    x_{norm} = K^{-1}x_{pixel} \approx [I|0]X \\
    x_{norm} = K'^{-1}x'_{pixel} \approx [R|t]X
\end{align}

We see that the world coordinate X is equivalent to the transpose of (x, 1). Then, if we perform the necessary matrix multiplication, we see that x' is equivalent to Rx + t. This shows us that there exists a linear relationship between x', Rx, and t, and thus all three lie on the same plane. The epipolar constraint can thus be written as the triple product of these, which is... 

\begin{align}
    x' \approx Rx + t \\ 
    x' \cdot [t \times (Rx)] = 0
\end{align}



\begin{align}
    x' \cdot [t \times (Rx)] = 0 \\ 
    x'^T[t]_{\times} Rx = 0
\end{align}

\begin{equation}
a \times b = 
    \begin{bmatrix}
        0 & -a_3 & a_2 \\
        a_3 & 0 & -a_1 \\
        -a_2 & a_1 & 0
    \end{bmatrix}
    \begin{pmatrix}
        b_1   \\
        b_2 \\
        b_3
    \end{pmatrix}
    = [a]_{\times}b
\end{equation}

As we know, the dot product of two vectors $a$ and $b$ is equivalent to $a$ transpose $b$, so therefore we can simplify the expression to the following:

\begin{align}
    x' \cdot [t \times (Rx)] = 0 \\ 
    x'^T[t]_{\times} Rx = 0 \\
    x'^TEx = 0 \\
    (x', y', 1)
    \begin{bmatrix}
        e_{11} & e_{12} & e_{13} \\
        e_{21} & e_{22} & e_{23} \\
        e_{31} & e_{32} & e_{33}
    \end{bmatrix}
    \begin{pmatrix}
        x   \\
        y \\
        1
    \end{pmatrix}
    = 0
\end{align}

### The Essential Matrix

We define the essential matrix to be $E = [t]_{\times}R$. **We use the essential matrix to relate corresponding points in stereo images.**

![](https://i.imgur.com/mTPRIJZ.png)

We can now transform our intuitive geometric reasoning about epipolar lines into more rigorous algebraic reasoning. Recall that a line can be expressed as $ax + by + c = 0$ or $l^T\mathbf{x}=0$ where $l=(a,b,c)^T$ and $\mathbf{x}=(x,y,1)^T$ in homogenous coordinates. Furthermore, recall the the essential matrix equation $x'^TEx=0$. If $l'=Ex$ we have $x'^Tl'=0$, meaning $x'$ must lie on $l'$. This means $Ex$ is the epipolar line associated $x$. The analogous case, $E^Tx'$ being the epipolar line associated with $x'$, is true as well due to the symmetry of the images.

In summary, the essential matrix $E$ has the following properties:

- $x'^TEx = 0$ where $E = [t]_{\times}R$
- $Ex$ is the epipolar line associated with $x (l'=Ex)$
- $E^Tx'$ is the epipolar line associated with $x'(l'=E^Tx')$
- Epipoles are the null spaces of $E: Ee = 0$ and $E^Tx'=0$ (note: all epipolar lines pass through epipoles)
- Rank 2 with 2 non-zero identical singular values
- 5 degrees of freedom (3 for $R$, 3 for $t$, -1 for scale)

### Uncalibrated Case 

**An overview ...**
- By substituting  $R$ & $t$ (Rotation/Translation information) into the epipolar constraint we can recover the essential matrix (see math above)
- We can decompose it, using matrix properties such as skew-symmetry and othonormal features, which allows SVD (Singular Value Decomposition)
- Since $x$ & $x'$ are 3D position of points that we dont have... we have an issue... BUT we have the projections of the scene points on the images. 
- We rewrite the epipolar constraint written with respect to the essential matrix in terms of the image coordinates
- The transformation of the essential matrix with the camera calibration matrices  yield the fundamental matrix which is written with respect to the image coordinates corresponding to the 3D point of interest and is implanted into our epipolar constraint. 
- Resolving the fundamental matrix from this constrain enables us to find the essential matrix and thus the Rotation and Translation matrices

**The math...**
![](https://i.imgur.com/MbW7MVM.png)

* We don't know the camera calibration matrices $K$ & $K'$
* Rewrite constraints in terms of *unknown normalized coordinates*

\begin{equation}
    x'^T_{norm} Ex_{norm} = 0
\end{equation}

where, 

\begin{align}
x_{norm} = K^{-1}x \\
x'_{norm} = K'^{-1}x'
\end{align}

As a result, $x'^T_{norm} Ex_{norm} = 0$ gives us;

\begin{equation}
x'^TFx = 0
\end{equation}

where we get $F$ (*the fundamental matrix*) discovered by Olivier Faugeras and Richard Hartley,

\begin{equation}
F=K'^{-T}EK^{-1}
\end{equation}
![](https://i.imgur.com/MeFqXD7.png)


* The fundamental matrix $F$ is a 3x3 matrix (note: everything is in homogenous coordinates) so **the fundamental matrix relates 2D homogenous coordinates in stereo images** where

\begin{equation}
    x'^TFx = 0 
\end{equation}

and can be expanded to the following representation:
\begin{align}
    (x', y', 1)
    \begin{bmatrix}
        f_{11} & f_{12} & f_{13} \\
        f_{21} & f_{22} & f_{23} \\
        f_{31} & f_{32} & f_{33}
    \end{bmatrix}
    \begin{pmatrix}
        x   \\
        y \\
        1
    \end{pmatrix}
    = 0
\end{align}
![](https://i.imgur.com/mTPRIJZ.png)

The fundamental matrix has the following properties:

- $Fx$ is the epipolar line assoced with $x(l' = Fx)$ which means $x'$ lies on the epipolar line
- $F^Tx'$ is the epipolar line associated with $x'(l=F^Tx')$
- Epipoles are the null spaces of $F$: $Fe=0$ and $F^Te'=0$ (note: the epipoles can be solved for via SVD)
- Rank 2 and singular
- 7 degrees of freedom (9 parameters, -1 for scale, -1 for singularity)

**Review of Essential and Fundamental Matrices**
- The *essential matrix* comes for the calibrated case (we know $K$ & $K'$) but the *fundamental matrix* comes about for the uncalibrated case when we do not know $K$ & $K'$
- Estimating the fundamental matrix is known as "weak calibration"
- Given the calibration matrix the essential matrix can be estimated
- The essential matrix provides the relative rotation and translation between the cameras
