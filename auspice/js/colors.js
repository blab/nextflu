var colors = ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"];
var regionColors = ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"]
var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]

var epitopeColorScale = d3.scale.linear().clamp([true])
	.domain(epiColorDomain)
	.range(colors);		

var nonepitopeColorScale = d3.scale.linear().clamp([true])
	.domain(nonEpiColorDomain)
	.range(colors);

var receptorBindingColorScale = d3.scale.linear().clamp([true])
	.domain(rbsColorDomain)
	.range(colors.filter( function(d,i){return i%2;}));

var lbiColorScale = d3.scale.linear()
	.domain([0.0, 0.02, 0.04, 0.07, 0.1, 0.2, 0.4, 0.7, 0.9, 1.0])
	.range(colors);

var dfreqColorScale = d3.scale.linear()
	.domain(dfreqColorDomain)
	.range(colors);

var HIColorScale_valid = d3.scale.linear()
	.domain(HIColorDomain)
	.range(colors);

var HIColorScale = function(c){
	if (c!='NaN'){
		return HIColorScale_valid(c);
	}else{
		return '#EEEEEE';
	}
}
HIColorScale.domain = HIColorScale_valid.domain;

var regionColorScale = d3.scale.ordinal()
	.domain(regions)
	.range(regionColors);

// "ep", "ne" and "rb" need no adjustments
function adjust_coloring_by_date() {
	if (colorBy == "lbi") {
		calcLBI(rootNode, nodes, false);
		nodes.forEach(function (d) {
			d.coloring = d.LBI;
		});
	}
	else if (colorBy == "dfreq") {
		calcDfreq(rootNode, freq_ii);
		nodes.forEach(function (d) {
			d.coloring = d.dfreq;
		});
	}
}


function colorByTrait() {
	
	colorBy = document.getElementById("coloring").value;
	console.log(colorBy);

	if (colorBy == "ep") {
		colorScale = epitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ep; });
	}
	else if (colorBy == "ne") {
		colorScale = nonepitopeColorScale;
		nodes.map(function(d) { d.coloring = d.ne; });
	}
	else if (colorBy == "rb") {
		colorScale = receptorBindingColorScale;
		nodes.map(function(d) { d.coloring = d.rb; });
	}
	else if (colorBy == "lbi") {
		colorScale = lbiColorScale;
		adjust_coloring_by_date();
	}
	else if (colorBy == "dfreq") {
		colorScale = dfreqColorScale;
		adjust_coloring_by_date();
	}
	else if (colorBy == "region") {
		colorScale = regionColorScale;
		nodes.map(function(d) { d.coloring = d.region; });
	}
	else if (colorBy == "cHI") {
		colorScale = HIColorScale;
		nodes.map(function(d) { d.coloring = d.cHI; });
	}

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
		
	d3.selectAll(".tip")
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
		
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();	 				
}

function tipStrokeColor(d) {
	var col = colorScale(d.coloring);
	return d3.rgb(col).toString();
}

function tipFillColor(d) {
	var col = colorScale(d.coloring);	;
	return d3.rgb(col).brighter([0.65]).toString();
}

function branchStrokeColor(d) {
	var col;
	if (colorBy == "region") {
		col = "#AAA";
	}
	else {
		col = colorScale(d.target.coloring);	
	}
	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
	return d3.rgb(modCol).toString();
}

function colorByGenotype() {
	var positions_string = document.getElementById("gt-color").value.split(',');
	var positions_list = []
	positions_string.map(function(d) {
		val = parseInt(d)-1;
		if (!isNaN(val)) {
			if (val < 551) {
				positions_list.push(val);
			}
		}
	});
	console.log(positions_list);
	if (positions_list.length > 0) {
		colorBy = "genotype";
		colorByGenotypePosition(positions_list);
	}
	else {
		d3.select("#coloring").each(colorByTrait);
	}
}

function colorByGenotypePosition (positions) {
	var gts = nodes.map(function (d) {
		var tmp = [];
		for (var i=0; i<positions.length; i++){
			var aa = cladeToSeq[d.clade];
			tmp[tmp.length] = (positions[i]+1)+aa[positions[i]];
		}
		d.coloring = tmp.join(" / "); 
		return d.coloring;});
	var unique_gts = d3.set(gts).values();
	var gt_counts = {};
	for (var i=0; i<unique_gts.length; i++){gt_counts[unique_gts[i]]=0;}
	gts.forEach(function (d) {gt_counts[d]+=1;});
	var filtered_gts = unique_gts.filter(function (d) {return gt_counts[d]>=10;});
	filtered_gts.sort(function (a,b){
		var res;
		if (gt_counts[a]>gt_counts[b]){ res=-1;}
		else if (gt_counts[a]<gt_counts[b]){ res=1;}
		else {res=0;}
		return res;});
	console.log("genotypes passed filtering:"+filtered_gts);
	colorScale = d3.scale.ordinal()
		.domain(filtered_gts)
		.range(genotypeColors);			
	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
	treeplot.selectAll(".tip")
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();
}

function colorByHIDistance(){
	correctVirus = document.getElementById("virus").checked;
	correctPotency = document.getElementById("serum").checked;
	predictedHI = document.getElementById("HIPrediction").checked;
	if (typeof(focusNode)=="undefined"){
		focusNode=rootNode;
	}
	treeplot.selectAll(".serum")
		.style("fill", function (d){if (d==focusNode) {return '#FF3300';} else {return '#555555';}})
		.style("font-size", function (d) {if (d==focusNode) {return "32px";} else {return "24px";}})
		.text(function (d) {if (d==focusNode) {return '\uf05b';} else {return '\uf10c';}});

	if (predictedHI){
		calcHIpred(focusNode, rootNode);
	}
	else{
		calcHImeasured(focusNode, rootNode);
	}
	console.log("Color by HI Distance from "+focusNode.strain);
	console.log("Using predictedHI: "+predictedHI);
	console.log("correcting for virus effect: "+correctVirus);
	console.log("correction for serum effect: "+correctPotency);

	colorScale = HIColorScale;
	nodes.map(function(d) { d.coloring = d.HI_dist;});

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
		
	d3.selectAll(".tip")
		.style("visibility", tipHIvalid)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
		
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();	 				
}


d3.select("#coloring")
	.style("cursor", "pointer")
	.on("change", colorByTrait);


d3.select("#gt-color")
	.on("keyup", colorByGenotype);



