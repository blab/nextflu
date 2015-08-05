Clazz.declarePackage ("J.renderspecial");
Clazz.load (["J.render.ShapeRenderer"], "J.renderspecial.PolyhedraRenderer", ["JU.P3", "JM.Atom", "JU.C"], function () {
c$ = Clazz.decorateAsClass (function () {
this.drawEdges = 0;
this.isAll = false;
this.frontOnly = false;
this.screens3f = null;
this.scrVib = null;
this.vibs = false;
Clazz.instantialize (this, arguments);
}, J.renderspecial, "PolyhedraRenderer", J.render.ShapeRenderer);
Clazz.overrideMethod (c$, "render", 
function () {
var polyhedra = this.shape;
var polyhedrons = polyhedra.polyhedrons;
this.drawEdges = polyhedra.drawEdges;
this.g3d.addRenderer (1073742182);
this.vibs = (this.ms.vibrations != null && this.tm.vibrationOn);
var needTranslucent = false;
for (var i = polyhedra.polyhedronCount; --i >= 0; ) if (polyhedrons[i].isValid && this.render1 (polyhedrons[i])) needTranslucent = true;

return needTranslucent;
});
Clazz.defineMethod (c$, "render1", 
 function (p) {
if (p.visibilityFlags == 0) return false;
var colixes = (this.shape).colixes;
var iAtom = p.centralAtom.i;
var colix = (colixes == null || iAtom >= colixes.length ? 0 : colixes[iAtom]);
colix = JU.C.getColixInherited (colix, p.centralAtom.colixAtom);
var needTranslucent = false;
if (JU.C.renderPass2 (colix)) {
needTranslucent = true;
} else if (!this.g3d.setC (colix)) {
return false;
}var vertices = p.vertices;
if (this.screens3f == null || this.screens3f.length < vertices.length) {
this.screens3f =  new Array (vertices.length);
for (var i = vertices.length; --i >= 0; ) this.screens3f[i] =  new JU.P3 ();

}var sc = this.screens3f;
var planes = p.planes;
for (var i = vertices.length; --i >= 0; ) {
var atom = (Clazz.instanceOf (vertices[i], JM.Atom) ? vertices[i] : null);
if (atom == null) {
this.tm.transformPtScrT3 (vertices[i], sc[i]);
} else if (atom.isVisible (this.myVisibilityFlag)) {
sc[i].set (atom.sX, atom.sY, atom.sZ);
} else if (this.vibs && atom.hasVibration ()) {
this.scrVib = this.tm.transformPtVib (atom, this.ms.vibrations[atom.i]);
sc[i].set (this.scrVib.x, this.scrVib.y, this.scrVib.z);
} else {
this.tm.transformPt3f (atom, sc[i]);
}}
this.isAll = (this.drawEdges == 1);
this.frontOnly = (this.drawEdges == 2);
if (!needTranslucent || this.g3d.setC (colix)) for (var i = planes.length; --i >= 0; ) {
var pl = planes[i];
this.fillFace (p.normixes[i], sc[pl.x], sc[pl.y], sc[pl.z]);
}
if (p.colixEdge != 0) colix = p.colixEdge;
if (this.g3d.setC (JU.C.getColixTranslucent3 (colix, false, 0))) for (var i = planes.length; --i >= 0; ) {
var pl = planes[i];
this.drawFace (p.normixes[i], sc[pl.x], sc[pl.y], sc[pl.z]);
}
return needTranslucent;
}, "J.shapespecial.Polyhedron");
Clazz.defineMethod (c$, "drawFace", 
 function (normix, a, b, c) {
if (this.isAll || this.frontOnly && this.vwr.gdata.isDirectedTowardsCamera (normix)) {
this.drawCylinderTriangle (a, b, c);
}}, "~N,JU.P3,JU.P3,JU.P3");
Clazz.defineMethod (c$, "drawCylinderTriangle", 
 function (a, b, c) {
var d = (this.g3d.isAntialiased () ? 6 : 3);
this.g3d.fillCylinderBits (3, d, a, b);
this.g3d.fillCylinderBits (3, d, b, c);
this.g3d.fillCylinderBits (3, d, a, c);
}, "JU.P3,JU.P3,JU.P3");
Clazz.defineMethod (c$, "fillFace", 
 function (normix, a, b, c) {
this.g3d.fillTriangleTwoSided (normix, a, b, c);
}, "~N,JU.P3,JU.P3,JU.P3");
});
