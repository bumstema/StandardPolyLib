(* ::Package:: *)

BeginPackage["StandardPolyLib`"]
ListOfPolygons={
(* Visualization Tools - Plots *)
PolyLib`StandardPolygon::usage,
PolyLib`SuperDisk::usage,
PolyLib`RadialFlower::usage,
PolyLib`StandardRectangle::usage,
PolyLib`Ellipse::usage,
PolyLib`HairyEllipse::usage,
PolyLib`HyperEllipse::usage,
PolyLib`SuperEllipse::usage,
PolyLib`Epicycloid::usage,
PolyLib`Hypocycloid::usage,
PolyLib`ReuleauxPolygon::usage,
PolyLib`Birkenstock::usage,
PolyLib`AlhambraCurve::usage,
PolyLib`AlhambraNail::usage,
PolyLib`StarPolygon::usage,
PolyLib`Heart::usage,
PolyLib`Butterfly::usage,
PolyLib`Pill::usage,
PolyLib`Leaf::usage,
PolyLib`RandomLeaf::usage,
PolyLib`KoshSnowflake::usage,
PolyLib`MinkowskiFractal::usage,
PolyLib`DragonFractal::usage
};

Render::usage=
"Render[poly_, colour___]: Quick display of polygon"
PolygonCentroid::usage=
"PolygonCentroid[polygon]: Gives the Center Of Mass for the Polygon"
inPolyQ2forPolyLib::usage=
"Returns 'True' or 'False' depending on if it is inside the polygon.\[IndentingNewLine]http://mathematica.stackexchange.com/questions/9405/how-to-check-if-a-2d-point-is-in-a-polygon"
OpenThePolygon::usage=
"OpenThePolygon[polygon_]:  If the first and last points are close, it removes the last point"
CloseThePolygon::usage=
"CloseThePolygon[polygon_]:  If the first and last points not close, it adds the last point"

MaxRadiusOfPolygon::usage="MaxRadiusOfPolygon[polygon]: Gives the maximum distance a point defining the polygon is from the origin {0,0}."
NormalizePolygon::usage="NormalizePolygon[poly_]: Centers the polygon at (0,0) and normalizes it to have the longest side have length 1.0."

StandardPolygon::usage= "StandardPolygon[ n_ ]:  Makes a Standard Polygon with 'n' Points (Circle as n -> \[Infinity])."
SuperDisk::usage="SuperDisk[p_,resolution_]: Uses  r = (\!\(\*SuperscriptBox[\(x\), \(p/2\)]\)+ \!\(\*SuperscriptBox[\(y\), \(p/2\)]\))  with resolution number of points."
RadialFlower::usage="RadialFlower[lobes_,fatness_,resolution_]: Uses  r = (fat + Cos[lobes]) on Circle with resolution number of points."
StandardRectangle::usage=
"StandardRectangle[b_]: Uses the equation of an ellipse to define the ratio (b) at 4 points =  \[Pi]/4 + n\[Pi]/2"
Ellipse::usage="Ellipse[b_,resolution_]: Uses  r = (b(\!\(\*SuperscriptBox[\(x\), \(2\)]\))+ \!\(\*SuperscriptBox[\(y\), \(2\)]\))  with resolution number of points. Uses x as the long axis."
HairyEllipse::usage="HairyEllipse[b_,hairLength_,resolution_]: Uses  r = (b(\!\(\*SuperscriptBox[\(x\), \(2\)]\))+ \!\(\*SuperscriptBox[\(y\), \(2\)]\)) but adds a hair with 'hairLength' radius to the ellipse at every vertex with resolution number of points."
HyperEllipse::usage=
"HyperEllispe[a_,p_,resolution_]: Builds a Hyper Ellipse, similar to the SuperDisk.   Uses  r = ( (a)\!\(\*SuperscriptBox[\(x\), \(p/2\)]\)+ \!\(\*SuperscriptBox[\(y\), \(p/2\)]\)) with x as the long axis"
SuperEllipse::usage=
"SuperEllipse[a_,b_,m_,n1_,n2_,n3_,resolution_]:  Creates a Polygon described in 'A GENERIC GEOMETRIC TRANSFORMATION THAT UNIFIES A WIDE RANGE OF NATURAL AND ABSTRACT SHAPES', American Journal of Botany 90(3): 333\[Dash]338. 2003. By JOHAN GIELIS"
Epicycloid::usage=
"Epicycloid[lobes_,resolution_]: Creates Epicycloid (Cloud-Shape) with specified number of lobes with resolution number of points."
Hypocycloid::usage=
"Hypocloid[lobes_,resolution_]: Creates Hypocycloid (Gear-Shape) with specified number of lobes with resolution number of points.  Note: Cannot have less than 3 Lobes."
ReuleauxPolygon::usage=
"ReuleauxPolygon[edges_, resolution_]: Builds the Reuleaux Polygon with edges number of sides and resolution number of points.  Will return a polygon with Mod[edges] points"
Birkenstock::usage=
"Birkenstock[p_, resolution_]: Builds Axe head shape with p width and Mod[4] resolution number of points."
AlhambraCurve::usage=
"AlhambraCurve[sides_,p_,resolution_]: Builds the Alhambra Curve with sides number of lobes, using p as the width (typically less than 1), and resolution Mod[2*sides] number of points."
AlhambraNail::usage=
"AlhambraNail[p_]: Builds the Alhambra Nail using p as the width (typically less than 1), and 9 points."
StarPolygon::usage=
"StarPolygon[points_, width___]:  Builds a star shape using different radius of a circle with 'points' number of points.  'width' can be specified (but not necessary with Default = 0.7 )to determine depth, R = 1+width at points, R=width at inner points"
Heart::usage=
"Heart[resolution_]:  Builds a Heart with resolution Points"
Butterfly::usage=
"Butterfly[resolution_]: Builds a butterfly with resolution Points"
Pill::usage="Pill[ nRadius_ ,  resolution_]:  Uses the Standard Polygon with an elongated middle section.  The space between the radius of the two end circles is nRadius."
Leaf::usage="Leaf[type_,resolution_]: Creates a leaf shape.  Only type=1 or type=2 are implemented leaves."
RandomLeaf::usage="RandomLeaf : Generates a leaf with random parameters with 200 points"
KoshSnowflake::usage="KoshSnowflake[ksidenum_,kiterations_,invert_]: Builds a Kosh Snowflake, where ksidenum is the number of sides, kiterations is the number of repeats, and invert is a boolean (True/False) that inverts the routine."
MinkowskiFractal::usage="MinkowskiFractal[msidenum_,miterations_]:= Returns the Minkowski fractal polygon, with msidenum being the symmetry number of sides, and miterations being the number of repeats."
DragonFractal::usage="DragonFractal[n_,inlev_,d_] := Builds a Dragon Type Polygon, with n being the number of legs, inlev being the number of iterations, and d being a squishing factor.  If[d>n, use n instead of d]."

Begin["PrivatePoly`"]
Render[poly_ , colour___]:=Graphics[{EdgeForm[Thick],If[colour===Null,Gray,colour],Polygon[poly]}];
PolygonCentroid[pts_?MatrixQ]:=With[{dif=Map[Det,Partition[pts,2,1]]},((ListConvolve[{1,1},#]&/@Transpose[pts]).dif)/(3 Total[dif])]
(* Returns "True" or "False" depending on if it is inside the polygon *)
(* http://mathematica.stackexchange.com/questions/9405/how-to-check-if-a-2d-point-is-in-a-polygon  *)
inPolyQ2forPolyLib = Compile[{{poly,_Real,2},{x,_Real},{y,_Real}},
Block[{Xi,Yi,Xip1,Yip1,u,v,w}, 
{Xi,Yi}=Transpose@poly;
Xip1=RotateLeft@Xi;
Yip1=RotateLeft@Yi;
u=UnitStep[y-Yi];
v=RotateLeft@u;
w=UnitStep[-((Xip1-Xi)(y-Yi)-(x-Xi)(Yip1-Yi))];
Total[(u (1-v)(1-w) - (1-u) v w)] !=0],CompilationTarget->"C"];
OpenThePolygon[poly_]:=Block[{},
If[Total@Abs[poly[[Length[poly] ]] -poly[[1]]  ] <= 10.^-3,Drop[poly,{Length[poly]}], poly]
];


CloseThePolygon[poly_]:=Block[{},
If[Total@Abs[poly[[Length[poly] ]] -poly[[1]]  ] >=  10.^-3,AppendTo[poly,poly[[1]] ]]
];

MaxRadiusOfPolygon[poly_]:=Block[{center,point,p,max},
center={0,0};
Max@Table[EuclideanDistance[center,poly[[p]]  ],{p,1,Length[poly]}]
];

NormalizePolygon[poly_]:=Block[{center,norm},
norm=poly;
Do[
center=PolygonCentroid[norm];
norm=Table[norm[[x]]-center,{x,Length[norm]}];
norm= N[norm/MaxRadiusOfPolygon[norm]];
,{x,1,12}];
norm
]
StandardPolygon[n_]:=N@Table[{Cos[x],Sin[x]},{x,0,2Pi,2Pi/n}];
SuperDisk[p_,resolution_]:=Block[{aa,ab,ia,ib,ic,id,int,r},
(*The radius as a function of the angle*)
r[\[Theta]_]:=1/Sqrt[(Sin[\[Theta]]^2)^p+(Cos[\[Theta]]^2)^p]^(1/p);
If[p<1,
(*Depending on the shape parameter, the density of points must be higher at points with higher curvature. Therefore we cannot use a linear angle dependence*)
(*ab and aa cover an angular range of \[Pi]/2 from -\[Pi]/4 to \[Pi]/4*)
ab=Table[\[Pi]/4 (t^(1/4)+1),{t,0,1-8/resolution,8/resolution}];
aa=Table[\[Pi]/4 (-t^(1/4)+1),{t,1,8/resolution,-8/resolution}];
ia=Join[aa,ab];
(*Increase the angles successively by \[Pi]/2*)
ib=ia+\[Pi]/2;
ic=ib+\[Pi]/2;
id=ic+\[Pi]/2;
int=Join[ia,ib,ic,id];,
ia=Table[2\[Pi] (t-1/2)^3+\[Pi]/4,{t,0,1-4/resolution,4/resolution}];
ib=ia+\[Pi]/2;
ic=ib+\[Pi]/2;
id=ic+\[Pi]/2;
int=Join[ia,ib,ic,id];
];
(*Generate the list of points*)
N[Table[{r[\[Theta]] Cos[\[Theta]],r[\[Theta]] Sin[\[Theta]]},{\[Theta],int}]]
];
RadialFlower[lobes_,fatness_,resolution_]:=Block[{rf},
rf= Table[ (fatness +Cos[lobes  \[Theta]]) { Cos[\[Theta]] , Sin[\[Theta]] } ,{\[Theta],0,2Pi,2Pi/resolution}];
N[rf/Max[rf]]
];
StandardRectangle[b_]:=Block[{rectangle},
rectangle=Table[{b Cos[\[Theta]],Sin[\[Theta]]},{\[Theta],\[Pi]/4,7\[Pi]/4, \[Pi]/2}];
N[rectangle/Max[rectangle]]
];
Ellipse[b_, resolution_]:=Block[{ellipse},
ellipse=Table[{b Cos[\[Theta]],Sin[\[Theta]]},{\[Theta],0,2Pi,2Pi/resolution}];
N[ellipse/Max[ellipse]]
];
HairyEllipse[b_,hairLength_, resolution_]:=Block[{hairy,n},
hairy=Table[{},{\[Theta],0,2Pi,2Pi/(3resolution)}];
n=1;
Do[
hairy[[n]]={b Cos[\[Theta]],Sin[\[Theta]]};
If[\[Theta]!= 2Pi,
hairy[[n+1]]=hairLength{b Cos[\[Theta]+10^-4],Sin[\[Theta]+10^-4]};hairy[[n+2]]={b Cos[\[Theta]+(2 10^-4)],Sin[\[Theta]+(2 10^-4)]}
];
n=n+3;,
{\[Theta],0,2Pi,2Pi/resolution}];
N[hairy/Max[hairy]]
];

(* Bjorn Implementation of Super Ellipse *)
HyperEllipse[a_, p_, resolution_]:=Block[{\[Theta],ab,aa,ia,ib,ic,id,int,poly,r},
r[\[Theta]_]=1/Sqrt[(Sin[\[Theta]]^2)^p+(Cos[\[Theta]]^2)^p]^(1/p);
If[p<1,
ab=Table[\[Pi]/4 (t^(1/4)+1),{t,0,1-8/resolution,8/resolution}];
aa=Table[\[Pi]/4 (-t^(1/4)+1),{t,1,8/resolution,-8/resolution}];
ia=Join[aa,ab];
ib=ia+\[Pi]/2;
ic=ib+\[Pi]/2;
id=ic+\[Pi]/2;
int=Join[ia,ib,ic,id];,
ia=Table[2\[Pi] (t-1/2)^3+\[Pi]/4,{t,0,1-4/resolution,4/resolution}];
ib=ia+\[Pi]/2;
ic=ib+\[Pi]/2;
id=ic+\[Pi]/2;
int=Join[ia,ib,ic,id];
];
poly=Table[{r[\[Theta]] Cos[\[Theta]],r[\[Theta]] Sin[\[Theta]]/a},{\[Theta],int}];
N[poly/Max[Abs@poly]]
];
SuperEllipse[a_,b_,m_, n1_, n2_,n3_, resolution_]:=Block[{\[Theta],r, poly,\[Phi],centerPoint},
r[\[Phi]_]:= N@(((Abs[Cos[\[Phi] m/4.]/a ]^n2)+ (Abs[Sin[\[Phi] m/4.]/b]^n3))^(-1./n1));
poly= Table[ r[\[Theta]]{Cos[\[Theta]],Sin[\[Theta]]},{\[Theta],0,2\[Pi], 2\[Pi]/resolution}];

NormalizePolygon[poly]
]

Epicycloid[lobes_,resolution_]:=Block[{epic},
(*Generate the list of points. For the formula, see wikipedia*)
epic=Table[{
(lobes+1)Cos[ \[Phi]]- Cos[(lobes+1) \[Phi]],
(lobes+1)Sin[ \[Phi]]- Sin[(lobes+1) \[Phi]]},
{\[Phi],0,2\[Pi], 2\[Pi]/resolution}];
N[epic/Max[Abs@epic]]
];
Hypocycloid[lobes_,resolution_]:=Block[{hypo},
hypo=Table[{
(lobes-1)Cos[ \[Phi]]+ Cos[(lobes-1) \[Phi]],
(lobes-1)Sin[ \[Phi]]- Sin[(lobes-1) \[Phi]]},
{\[Phi],0,2\[Pi], 2\[Pi]/resolution}];
N[hypo/Max[Abs@hypo]]
];
ReuleauxPolygon[edges_, resolution_]:=Block[{corners,ppe,\[Alpha],r,points,n,m},
(*The location of the corners*)
corners=Table[{Cos[2m \[Pi]/edges],Sin[2m \[Pi]/edges]},{m,0,edges-1}];
(*The number of points per edge (ppe)*)
ppe=Ceiling[resolution/edges];
(*The half opening angle of the circle segments*)
\[Alpha]=ArcTan[Sin[\[Pi]/edges]/(Cos[\[Pi]/edges]+1)];
(*The radius of the circle segments*)
r=N@Sqrt[(Cos[\[Pi]/edges]+1)^2+Sin[\[Pi]/edges]^2];
(*The list of vertices*)
points={};
For[m=1,m<=Length[corners],m++,
(*Starting from each corner, a circle segment is represented by line segments*)
For[n=0,n<ppe,n++,
AppendTo[points,
corners[[m]]+
{r Cos[2(m-1) \[Pi]/edges+\[Pi]-\[Alpha]+2n \[Alpha]/ppe],
r Sin[2(m-1)\[Pi]/edges+\[Pi]-\[Alpha]+2n \[Alpha]/ppe]
}]
];
];
N[points/Max[points]]
];
Birkenstock[p_, resolution_]:=Block[{\[Alpha],r,resa,p0,arc1,a0,arc2,m,curve1,curve2,poly,centerPoint},
(*Generate the list of points, with p=b/a (a:=1) being the parameter*)
(*Determine the half angle of one arc*)
\[Alpha]=ArcTan[1/(2p)];
(*Determine the radius of one arc*)
r=Sqrt[1/4+p^2];
(*The resolution of one arc*)
resa=(resolution/4);
(*The center of the right arc*)
p0={1/2-p,0};
(*Determine first arc*)
arc1=Table[p0+{r*Cos[-\[Alpha]+\[Phi]],r*Sin[-\[Alpha]+\[Phi]]},{\[Phi],0,2\[Alpha],(2\[Alpha])/resa}];
a0={1/2-p,0};
(*Translate first arc*)
arc2=MapIndexed[Subtract[#1,a0]&,arc1[[2;;Length[arc1]-1]]];
m=RotationMatrix[-(\[Pi]/2)];
(*Copy first artc partially and rotate by \[Pi] and reverse the order of the points*)
arc2=Reverse[Dot[m,#1]&/@arc2];
a0={0,1/2+p};
arc2=MapIndexed[Plus[#1,a0]&,arc2];
(*Joins both arcs and one "half" is finished*)
curve1=Join[arc1,arc2];
(*Rotate this half by 180 degrees*)
m=RotationMatrix[\[Pi]];
curve2=Dot[m,#1]&/@curve1;
(*Join the curves and voila*)
poly=Join[curve1,curve2];

centerPoint=PolygonCentroid[poly];
poly=poly[[#]]-centerPoint&/@Range[Length[poly]];

N[poly/Max[Abs@poly]]
];
AlhambraCurve[sides_,p_,resolution_]:=Block[{\[Alpha],r,resa,p0,arc1,a0,m,arc2,curve,ri,curve1,rm,curves,poly},
(*Generate the list of points, with p=b/a (a:=1) being the parameter*)
(*Determine the half angle of one arc*)
\[Alpha]=ArcTan[1/(4p) ];
(*Determine the radius of one arc*)
r=N@Sqrt[1/16+p^2];
(*The resolution of one arc*)
resa=resolution/(2 sides);
(*The number of sides of the Alhambra curve*)
p0={1/4,-p};
(*Determine first arc*)
arc1=Table[p0+{r*Cos[\[Pi]/2-\[Alpha]+\[Phi]],r*Sin[\[Pi]/2-\[Alpha]+\[Phi]]},{\[Phi],0,2\[Alpha],(2\[Alpha])/resa}];
a0=arc1[[1]];
(*Translate first arc*)
arc1=MapIndexed[Subtract[#1,a0]&,arc1];
m=RotationMatrix[\[Pi]];
(*Copy first artc partially and rotate by \[Pi] and reverse the order of the points*)
arc2=arc1[[2;;Length[arc1]-1]];
arc2=Reverse[Dot[m,#1]&/@arc2];
(*Joins both arcs and the side is finished*)
curve=Join[arc2,arc1];
(*Create sides sides and connect them*)
(*Determine the inner radius of the n-edge*)
ri=1/(2Tan[\[Pi]/sides]);
(*Move curve up by ri*)
curve1=MapIndexed[Plus[#1,{0,ri}]&,curve];
(*Create the rotation matrices for each side*)
rm=Table[RotationMatrix[n (2\[Pi])/sides],{n,1,sides-1}];
(*Rotate curve1 by each rotation matrix and create a list of sides*)
curves=Table[Dot[rm[[n]],#1]&/@curve1,{n,1,Length[rm]}];
(*Join the sides*)
poly=Join[curve1,Flatten[curves,1]];
N[poly/Max[Abs@poly]]
];

AlhambraNail[p_]:=Block[{\[Theta],points,m},
\[Theta]=ArcTan[(1-p)/(1+p)];
points={{0,0},{0,1},{p,1},{p,1+p},{-(1-p),1+p}};
(*Rotation Matrix -\[Theta] and rotation if the RHS*)
m=RotationMatrix[-\[Theta]];
points=Dot[m,#1]&/@points;
(*Add the points missing now on the LHS, be careful with the order of the points*)
points=Join[points,{{-points[[4,1]],points[[4,2]]},{-points[[3,1]],points[[3,2]]},{-points[[2,1]],points[[2,2]]}}];
AppendTo[points,points[[1]]];
N[points/Max[Abs@points]]
];
StarPolygon[points_, width___]:=Block[{star,w},
If[width===Null,w=0.7,w=width];
star=Table[  (Mod[\[Theta],2\[Pi]/(points)]+w){ Cos[\[Theta]],Sin[\[Theta]]},{\[Theta],0,2\[Pi], \[Pi]/points}];
N[star/Max[star]]
];
Heart[ resolution_]:=Block[{x,y,heart,centerPoint},

x[t_]:=16Sin[t]^3;
y[t_]:=13Cos[t]-5Cos[2t]-2Cos[3t]-Cos[4t] ;

heart=N@Table[{x[t],y[t]},{t,0,2Pi,2Pi/resolution}];

NormalizePolygon[heart]

];
Butterfly[ resolution_]:=Block[{x,y,butterfly,centerPoint, outerpoints},

x[u_]:=Cos[u ] (Exp[Cos[u]]- 2Cos[4u]-Sin[(u/12)]^(5.));
y[u_]:=Sin[u ] (Exp[Cos[u]]- 2 Cos[4 u]-Sin[(u/12)]^(5.));

butterfly=Table[{x[\[Theta]],y[\[Theta]]},{\[Theta],0,2Pi,2Pi/resolution}];

outerpoints= Table[
If[
inPolyQ2forPolyLib[Drop[butterfly,{x}],  butterfly[[x,1]],  butterfly[[x,2]] ] == False,
butterfly[[x]],
{}]
,{x,1,Length[butterfly]}];


outerpoints=DeleteCases[outerpoints,{}];

centerPoint=PolygonCentroid[outerpoints];
butterfly=outerpoints[[#]]-centerPoint&/@Range[Length[outerpoints]];

N[butterfly/Max[butterfly]]
];
Pill[ nRadius_ ,  resolution_]:= Block[{pill},
pill=Join[
Table[{Cos[x],Sin[x]}+{nRadius/2, 0},{x,0,Pi/2,2Pi/(resolution)}],
Table[{Cos[x],Sin[x]}-{nRadius/2, 0},{x,Pi/2,3Pi/2,2Pi/( resolution)}],
Table[{Cos[x],Sin[x]}+{nRadius/2,0},{x,3Pi/2,2Pi,2Pi/(resolution)}]
];

N[pill/Max[pill]]
];
Leaf[type_,resolution_]:=Block[{leaf,center},
If[type==1,
leaf=Table[(1+Sin[t])(1+0.3Cos[8t])(1+0.1Cos[24t]){Cos[t],Sin[t]},{t,0,2Pi,2Pi/resolution}];,
leaf=Table[(100./(100.+(t-Pi/2)^8))((2-Sin[7t]-Cos[30t]/2)){Cos[t],Sin[t]},{t,-Pi/2,3/2Pi, 2 Pi/resolution}];
];

NormalizePolygon[leaf]

];

RandomLeaf:=Block[{ randA,randB,randC,randD,randomLeaf,center},
 randA=RandomReal[];
randB=RandomReal[];
randC= RandomInteger[10];
randD= RandomInteger[24];
randomLeaf=Table[(1+Sin[t])(1+randA Cos[randC t])(1+randB Cos[randD t]){Cos[t],Sin[t]},{t,0,2Pi,Pi/100}];

NormalizePolygon[randomLeaf]
];
(*http://demonstrations.wolfram.com/FractalizingPolygons/*)
KoshSnowflake[ksidenum_,kiterations_,invert_]:=
Block[{points,tpoints,kosh},
points=Reverse[N[Table[{Cos[t*2Pi/#+Pi/2],Sin[t*2Pi/#+Pi/2]},{t,0,#}]]]&;
tpoints=points[ksidenum];

kosh=Nest[Flatten[Partition[#,2,1]/.{{a_,b_},{c_,d_}}:> {{a,b},{a+(c-a)/3,b+(d-b)/3},{a+(c-a)/2-(d-b)/(2*Sqrt[3]),b+(d-b)/2+(c-a)/(2*Sqrt[3])},{a+2*(c-a)/3,b+2*(d-b)/3},{c,d}},1]&,If[invert,Reverse,Identity][points[ksidenum]],kiterations];

NormalizePolygon[kosh]

];
(*http://demonstrations.wolfram.com/FractalizingPolygons/*)
MinkowskiFractal[msidenum_,miterations_]:=
Block[{points,tpoints,minkowski},
points=Reverse[N[Table[{Cos[t*2Pi/#+Pi/2],Sin[t*2Pi/#+Pi/2]},{t,0,#}]]]&;
tpoints=points[msidenum];

minkowski=Nest[Flatten[Partition[#,2,1]/.{{a_,b_},{c_,d_}}:>{{a,b},{a+(c-a)/4,b+(d-b)/4},{a+(c-a)/4-(d-b)/4,b+(d-b)/4+(c-a)/4},{a+2(c-a)/4-(d-b)/4,b+2(d-b)/4+(c-a)/4},{a+2(c-a)/4,b+2(d-b)/4},{a+2(c-a)/4+(d-b)/4,b+2(d-b)/4-(c-a)/4},{a+3(c-a)/4+(d-b)/4,b+3(d-b)/4-(c-a)/4},{a+3(c-a)/4,b+3(d-b)/4},{c,d}},1]&,points[msidenum],miterations];

NormalizePolygon[minkowski]
];
(*http://demonstrations.wolfram.com/TheBoundaryOfPeriodicIteratedFunctionSystems/*)
DragonFractal[n_,inlev_,d_]:=Block[{gr,dodaj,i1,i2,il,nr,mls,mx,mkline,pos,len,a,ln,trl,tr,ac,b,l,mr,dir,r,nl,DD,z,one, zl,zmd,uz,ie,dat,do,od,po,ko,ki,jj,j,nk,npo,fin,lam,v, w,i,maxlv, lev,lv,inln,lntb,nln,mne,ml,hor},
lev=inlev;
lv=lev;
If[d> n,DD=-n,DD=-d];

lam=2.1304;fin={{1,1,1,1,0,0},{0,0,0,0,1,0},{0,1,0,1,0,0},{1,1,0,0,0,1},{0,0,1,0,0,1},{0,1,1,0,0,0}};v={0.1303954347672788`,0.27779383894270526`,0.18362145258003185`,0.14739840417542643`,0.1303954347672788`,0.1303954347672788`};w={{2,0,2,0,0,0},{2,4,2,4,0,0},{6,8,2,4,4,4},{10,16,14,8,8,8},{18,40,26,24,16,16},{42,84,50,44,40,40},{86,176,122,92,84,84},{178,384,254,208,176,176},{386,816,530,432,384,384},{818,1732,1154,916,816,816},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};

gr={};
dodaj[a_,b_,c_]:={
i1=6*a+b+1;i2=6*Mod[a+3,6]+Mod[b+3,6]+1;If[a>2,il=i2,il=i1];If[(nr[[il]]==0)&&(c>0),uz[[++ie]]=il-1;nr[[i1]]=nr[[i2]]=ie;mls=Append[mls,{}]];
mx[[i,nr[[il]]]]+=c;mls[[i]]=Join[mls[[i]],Table[nr[[il]],{i,c}]]};
mkline[lv_,ps_,ca_]:={pos=ps;len=Length[l];If[ca,a=Mod[Quotient[l[[1]],6]+3,6]];ln=Table[pos,{k,len+1}];trl=tr[[lv]];
Do[ac=l[[k]];b=Mod[ac,6];If[a!=Quotient[ac,6],a=b,a=Mod[b+3,6]];ln[[k+1]]=(pos+=mr*trl[[a+1]]),{k,len}];};


z=DD/2+I*Sqrt[n-(DD/2)^2]+0.;r=N[Sqrt[n]];mr=1;
tr=Table[one={Re[zl=z^-lv],Im[zl]};zmd={Re[zl*(z-DD)],Im[zl*(z-DD)]};{one,one-zmd,-zmd,-one,zmd-one,zmd},{lv,1,20}];
mx=Table[0,{i,1,18},{j,1,18}];uz=Table[0,{i,1,18}];nr=Table[0,{i,1,36}];mls={{}};uz[[1]]=0;nr[[1]]=nr[[22]]=1;ie=1;
dat={{n+DD,0,4},{n-1,n+DD,5},{n-1,n-1,0},{-DD-1,n-1,1},{0,-DD-1,2},{0,0,3}};
For[i=1,i<=ie,++i,
do=Mod[a=uz[[i]],6];od=(a-do)/6;
{po,ko,ki}=dat[[od+1]];
If[do<=od,do+=6];
For[j=od+1,j<=do,j++,jj=Mod[j,6];npo=dat[[jj+1,2]];If[po!=npo,If[npo<po,nk=3,nk=0];dodaj[ki,nk,1];dodaj[0,3,Abs[npo-po]-1];ki=3-nk;po=npo]];dodaj[ki,dat[[jj+1,3]],1]];
fin=Table[mx[[i,j]],{i,1,ie},{j,1,ie}];{{lam},{v}}=N[Eigensystem[Transpose[fin],1]];v/=Sum[v[[i]],{i,1,ie}];
w=Table[0,{j,1,30},{i,1,ie}];w[[1,nr[[4]]]]=2*(n-2);w[[1,nr[[1]]]]=2;
i=1;While[Sum[w[[i,k]],{k,1,ie}]<5000,w[[i+1]]=w[[i]].fin;i++];maxlv=i-1;lev=Min[maxlv,lev];
inln=Table[nr[[4]],{i,n-1}];inln[[n-1]]=nr[[1]];lntb=Table[{},{i,1,maxlv}];
Do[lntb[[i]]=Table[uz[[inln[[j]]]],{j,1,Length[inln]}];
If[i<maxlv,nln={};Do[nln=Join[nln,mls[[inln[[j]]]]],{j,1,Length[inln]}];inln=nln],{i,1,maxlv}];
mne=Sum[w[[lev,k]],{k,1,ie}];
Do[ac=uz[[i]];pos={0,-2*i};hor={1,0};If[n<=3,mr=r^lv,mr=r^(lv-1)];ln={pos,pos+mr*tr[[lv,1+Mod[ac,6]]]};
ml=mls[[i]];l=Table[uz[[ml[[k]]]],{k,Length[ml]}];a=Mod[Quotient[ac,6]+1,6];
mkline[lv+1,pos+2.5*hor,False];
mr=1;
,{i,1,ie}];
l=Join[lntb[[lev]],lntb[[lev]]]; mkline[lev,{0,0},True];

Return[NormalizePolygon[ln]];
];


End[]
EndPackage[]
