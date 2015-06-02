var regions = [
//generic country codes
"GUI","GIN","SLE","LIB","LBR",
// Guinean prefectures
"Boffa", "Boke", "Fria", "Gaoual", "Koundara", "Conarky", "Dabola",
"Dinguiraye", "Faranah", "Kissidougou", "Kankan", "Kerouane",
"Kouroussa", "Mandiana", "Siguiri", "Coyah", "Dubreka", "Forecariah",
"Kindia", "Telimele", "Koubia", "Labe", "Lelouma", "Mali", "Tougue",
"Dalaba", "Mamou", "Pita", "Beyla", "Gueckedou", "Lola", "Macenta",
"Nzerekore", "Yamou",
// Sierra Leonean districts
"Kailahun", "Kenema", "Kono", "Bombali", "Kambia", "Koinadugu",
"PortLoko", "Tonkolili", "Bo", "Bonthe", "Moyamba", "Pujehun",
"WesternRural", "WesternUrban",
// Liberian counties
"Nimba", "RiverCess", "RiverGee", "Sinoe", "Bomi", "Bong",
"Gbapolu", "GrandCapeMount", "GrandBassa", "GrandGedeh",
"GrandKru", "Lofa", "Margibi", "Maryland", "Montserrado",
];
var restrictTo = "all";
var restrictToLab = "all";

var cladeToSeq = {}

var globalDate = new Date();

var nodes, tips, rootNode, links, vaccines;

var nDisplayTips, displayRoot;
if (document.getElementById("gtspec") != null){
    var freqdefault = document.getElementById("gtspec").value;
}else{
    var freqdefault ='';
}

function treePlotHeight(width) {
	return 400 + 0.35*width;
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
var epiColorDomain = genericDomain;
var nonEpiColorDomain = genericDomain;
var rbsColorDomain = genericDomain;
var dateColorDomain = genericDomain;
var dfreqColorDomain = genericDomain.map(function(d){return Math.round(100*(-0.18+d*0.36))/100;});
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
