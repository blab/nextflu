var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea", "NorthAmerica", "SoutheastAsia", "SouthAsia", "China"]

var cladeToSeq = {}

if (typeof globalDate == 'undefined') {
    var globalDate = new Date();
}

var nodes, tips, rootNode, links, vaccines, sera;

var nDisplayTips, displayRoot;
if (document.getElementById("gtspec") != null){
    var freqdefault = document.getElementById("gtspec").value;
}else{
    var freqdefault ='';
}

function treePlotHeight(width) {
	return 400 + 0.30*width;
}
var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
var treeWidth = containerWidth;
var treeHeight = treePlotHeight(treeWidth);
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
var serumSymbol = '\uf0fe';
var epiColorDomain = genericDomain;
var nonEpiColorDomain = genericDomain;
var rbsColorDomain = genericDomain;
var dateColorDomain = genericDomain;
var HIColorDomain = genericDomain.map(function(d){return Math.round(100*(d*3.6))/100;});
var dfreqColorDomain = genericDomain.map(function(d){return Math.round(100*(0.2+d*1.8))/100;});
var fitnessColorDomain = genericDomain.map(function(d){return Math.round(100*((d-0.5)*6.0))/100;});
var time_step;


d3.json(path + file_prefix + "meta.json", function(error, json) {
    if (error) return console.warn(error);
    d3.select("#updated").text(json['updated']);
    commit_id = json['commit'];
    short_id = commit_id.substring(0, 6);
    d3.select("#commit")
        .append("a")
        .attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
        .text(short_id);

});
