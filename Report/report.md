


# 1. Table detection

In order to take detect the table, we have to take advantage of its features. The table edges are long straight lines. It is also blue and stands out from its environment. So our approach consists of finding the lines present in the input image and then filter them according to different criteria. After that, the lines are grouped depending on which side of the table they belong to. From this group of lines, the four sides of the table is calculated

# 1.1 Lines detection

Using Hough peaks and Hough lines, the algorithm start by detecting the lines that are present in the input image. We consider as a valid lines the ones that are on the edge of the table. The lines are therefore filtered on the following criteria:

**Orientation**

Given the point of view of the camera, we don't expect to have an edge of the table showing vertically on the image. Hough lines algorithm uses a polar coordinates system, therefore we filter the lines directly from $\theta$ of the polar coordinates.

**Color separation**

A valid line will have from one side the blue table color and then from the other side any other color (floor green, wall white...). To decide which lines are valid, two sample pixels are taken from each side of the line. The midpoint of the line along with the orthogonal line going through the midpoint are calculated. The sample pixels are equally distant points from the centroid on the orthogonal lines. The HSV values of these points are considered as it is a more robust way to deduce the color.

Given two points $p_{1} = (h_{1}, s_{1}, v_{1})$ and $p_{2} = (h_{2}, s_{2}, v_{2})$, a color is considered blue if

$blue_{1} = (h_{1} > BLUEMIN) \wedge (h_{1} < BLUEMAX) \wedge (h_{1} > VLIGHTMIN)$

where the thresholds (BLUEMIN, BLUEMAX and VLIGHTMIN) are defined from the hue values of the blue color. A minimum value of brightness is also imposed as blacks may be considered as blue.

A line is considered valid if

$line = blue_{1} \veebar blue_{2}$  

## 1.2 Lines grouping

After detecting the lines that are on the edge of the table, we need to extract the four lines that outline the shape of the table. We start by grouping the lines into four groups. Using the k-means clustering algorithm, the lines are separated into four groups. The algorithm uses squared euclidean distance.

Given that each line is represented by $ax + b$, we calculated the weighted average of the lines. We gave more weight to the longest lines as they are more accurate.

The fact the we get many lines for each of the edges of the table is due to the radial distortion induced by the wide lens camera.

Calculating the corner points is then straight forward, as they are the intersection of the edge lines.

# 2. Camera Calibration and Position

## 2.1 Internal Parameters K

Knowing that the image of the absolute conic $\omega$ depends only on the internal parameters $K$ of the camera matrix $P$, we can find $K$ by calculating the image of the absolute conic $\omega$.

In addition to that, we know that two vanishing points $v_{1}$ and $v_{2}$ that correspond to perpendicular directions satisfy the following property:

$v_{1}^T\omega v_{2} = 0$

To generate enough constraints on $\omega$, we need to find multiple vanishing points corresponding to perpendicular directions. After detecting lines using Hough Lines, we choose the perpendicular lines and calculate their intersection with the vanishing line.

The vanishing line is found by: intersecting parallel lines to find vanishing points and the fitting a line going through these vanishing points found before.

Using the Choleski factorization, we compute the matrix $K$ from $\omega = (KK^T)^-1$

## 2.2 Cameras position

For each one of the cameras we start finding its localisation relatively to the table. We set the world reference frame equal to the camera reference frame.

We consider a table frame $O$, having its origin on the bottom left corner of the table with the $x$ axis along the short edge of the table and the $y$ axis along the long one. Any point on the table will be expressed as $X^{O} = \begin{pmatrix}
i^{w}_{O} & j^{w}_{O} & t^{w}_{O}
\end{pmatrix} = K^{-1}H^{T}$ in $O$. We can express the table reference frame according to the world reference frame as:

$X^{w}= \bigl(\begin{smallmatrix}
i^{w}_{O} & j^{w}_{O}  & k^{w}_{O} & t^{w}_{O} \\
0 & 0 & 0 & 1  
\end{smallmatrix}\bigr) X^{O}  =  \bigl(\begin{smallmatrix}
i^{w}_{O} & j^{w}_{O} & t^{w}_{O} \\
0 & 0 & 1  
\end{smallmatrix}\bigr) \begin{pmatrix}
x\\
y\\
w\\
\end{pmatrix}$

Using the homography that maps the table to real world coordinates and the internal parameters K we find the camera rotation and translation based on a pinhole model camera, we have:

$X'= \bigl(\begin{smallmatrix}
 & 0 \\
K &  0\\
 &  0
\end{smallmatrix}\bigr)\bigl(\begin{smallmatrix}
i^{w}_{O} & j^{w}_{O} & t^{w}_{O} \\
0 & 0 & 1  
\end{smallmatrix}\bigr) \begin{pmatrix}
x\\
y\\
w\\
\end{pmatrix} = K \bigl(\begin{smallmatrix}
i^{w}_{O} & j^{w}_{O} & t^{w}_{O}
\end{smallmatrix}\bigr)
 \begin{pmatrix}
x\\
y\\
w\\
\end{pmatrix}$

The homography mapping the table plan to the image coordinates has already been calculated and it satisfies:

$X' = H \begin{pmatrix}
x & y & w
\end{pmatrix} ^ {T}$ and therefore $\begin{pmatrix}
i^{w}_{O} & j^{w}_{O} & t^{w}_{O}
\end{pmatrix} = K^{-1}H$

$k^{w}_{O}$ can be found by calculating the cross product of $i^{w}_{O}$ and $j^{w}_{O}$. So the orientation of $O$ in the camera frame is $R^{w}_{O} = \begin{pmatrix} i^{w}_{O} & j^{w}_{O} & k^{w}_{O} \end{pmatrix}$

Using the same method, we can find $O'$, the second camera frame, and its localisation in the camera frame. And since we found the localisation of the first camera to the table and the table according to the second camera (inverse of the localisation of the second camera to the table), we can find the localisation of the cameras according to each other.

The 3D position of the ball in the table frame can be found using the two camera matrices.

## Epipolar geometry

The pinhole model used as a basis for the previous analysis is a 3D to 2D conversion. When doing the conversion we loose an important information which is the depth of the image. Therefore to find the depth, we use two (or more) cameras. The intrinsic projective geometry of two views is called epipolar geometry. It depends only on the internal parameters of the camera and their relative pose while being independent from the scene structure.

The main component of epipolar geometry is the $3\times 3$ matrix $F$ called the fundamental matrix and satisfies the following property:
$xFx^{'}=0$

Where $x$ and $x^{'}$ are the image coordinates of the real world point $X$ in the first and second view.

The points $x$, $x^{'}$, the camera centers ($c$, $c^{'}), and $X$ are coplanar and constitue the plane $\pi$ called the epipolar plane. The epipolar plane meets the 
