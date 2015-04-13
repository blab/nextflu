var restrictTo = "all";

var cladeToSeq = {}

var globalDate = new Date();

var nodes, tips, internals, rootNode, links, vaccines;

function treePlotHeight(width) {
	return 400 + 0.35*width;
}
var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
var treeWidth = containerWidth, treeHeight;
treeHeight = treePlotHeight(treeWidth);
var tree = d3.layout.tree()
	.size([treeHeight, treeWidth]);


var treeplot = d3.select("#treeplot")
	.attr("width", treeWidth)
	.attr("height", treeHeight);

var legend = d3.select("#legend")
	.attr("width", 280)
	.attr("height", 100);

var colorBy = document.getElementById("coloring").value;
var colorScale;

var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea", "NorthAmerica", "SoutheastAsia", "SouthAsia", "China"]

var pivots, dt;
