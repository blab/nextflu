var genericDomain = [0,0.111,0.222,0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 1.0];
var path = '/data/';
var tip_labels = true;

var countries = [];
var divisions = [];

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
var serumSymbol = '\uf013';
var epiColorDomain = genericDomain;
var nonEpiColorDomain = genericDomain;
var ageColorDomain = genericDomain.map(function(d){return Math.round(d*100);});
var ageScoreColorDomain = genericDomain.map(function(d){return Math.round(d*500)/100;});
var genderColorDomain = genericDomain.map(function(d){return Math.round((d-0.5)*200)/100;});
var rbsColorDomain = genericDomain;
var dateColorDomain = genericDomain;
var HIColorDomain = genericDomain.map(function(d){return Math.round(100*(d*3.6))/100;});
var dfreqColorDomain = genericDomain.map(function(d){return Math.round(100*(0.7+d*0.6))/100;});
var fitnessColorDomain = genericDomain.map(function(d){return Math.round(100*((d-0.5)*16.0))/100;});
var time_step;

d3.json(path + file_prefix + "meta.json", function(error, json) {
    if (error) return console.warn(error);
    update_date = json['updated'];
    d3.select("#updated")
      .append("span")
      .html("updated " + update_date);
    commit_id = json['commit'];
    short_id = commit_id.substring(0, 6);
    if (commit_id !== "unknown") {
      d3.select("#commit")
        .append("span")
        .html("and processed with commit ")
        .append("a")
        .attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
        .text(short_id);
    }

});

String.prototype.toTitleCase = function() {
	return this.replace(/\w+/g, function(txt){return txt.charAt(0).toUpperCase() + txt.substr(1).toLowerCase();});
}
