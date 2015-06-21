// 2 color	["#5097BA", "#DF4327"]

var colors = [
	[], 
	["#8EBC66"], 
	["#5EA9A1", "#E0A23A"],
	["#5CA7A4", "#A4BE56", "#E39B39"],
	["#5AA5A7", "#8ABB6A", "#C3BA46", "#E69237"],
	["#59A3AA", "#7CB879", "#AABD52", "#D4B13F", "#E68634"],
	["#57A1AD", "#73B583", "#98BD5E", "#C0BA47", "#DEA73C", "#E67B32"],
	["#559EB1", "#6DB28B", "#8CBB68", "#AFBD4F", "#CFB541", "#E39B39", "#E56E30"],
	["#539CB4", "#68AF92", "#83BA71", "#A2BE57", "#C1BA47", "#DAAC3D", "#E68E36", "#E35F2D"],
	["#519AB7", "#64AC99", "#7BB879", "#97BD5F", "#B5BD4C", "#CFB541", "#E0A13A", "#E67F33", "#E1512A"],
	["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"],
	["#4E95BD", "#5DA8A3", "#70B487", "#87BB6D", "#A0BE58", "#BABC4A", "#D0B440", "#E0A23A", "#E68634", "#E35F2D", "#DD3425"],
	["#4D92BF", "#5AA5A8", "#6BB18D", "#80B974", "#98BD5E", "#B1BD4E", "#C8B944", "#DAAC3D", "#E59738", "#E67732", "#E14F2A", "#DB2522"]
];
var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]

var epitopeColorScale = d3.scale.linear().clamp([true])
	.domain(epiColorDomain)
	.range(colors[10]);		

var nonepitopeColorScale = d3.scale.linear().clamp([true])
	.domain(nonEpiColorDomain)
	.range(colors[10]);

var receptorBindingColorScale = d3.scale.linear().clamp([true])
	.domain(rbsColorDomain)
	.range(colors[4]);

var lbiColorScale = d3.scale.linear()
	.domain([0.0, 0.02, 0.04, 0.07, 0.1, 0.2, 0.4, 0.7, 0.9, 1.0])
	.range(colors[10]);

var dfreqColorScale = d3.scale.linear()
	.domain(dfreqColorDomain)
	.range(colors[10]);

var HIColorScale_valid = d3.scale.linear()
	.domain(HIColorDomain)
	.range(colors[10]);

var cHIColorScale = d3.scale.linear()
	.domain(HIColorDomain)
	.range(colors[10]);


var HIColorScale = function(c){
	if (c!='NaN'){
		return HIColorScale_valid(c);
	}else{
		return '#EEEEEE';
	}
}
HIColorScale.domain = HIColorScale_valid.domain;

var regionColorScale = d3.scale.ordinal()
	.domain(regions.map(function(d){return d[0];}))
	.range(regions.map(function(d){return d[1];}));

var dateColorScale = d3.scale.linear().clamp([true])
	.domain(dateColorDomain)
	.range(colors[10]);

var regionColorScale = d3.scale.ordinal()
	.domain(regions.map(function(d){return d[0];}))
	.range(regions.map(function(d){return d[1];}));

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
	else if (colorBy == "date") {
		nodes.forEach(function (d) {
			d.coloring = d.num_date;
		});
	}	
}

function stateAtPosition(clade, gene, pos){
	if (typeof cladeToSeq[clade][gene][pos] == "undefined"){
		return cladeToSeq["root"][gene][pos];
	}else{
		return cladeToSeq[clade][gene][pos];
	}
}

function colorByTrait() {
	
	colorBy = document.getElementById("coloring").value;
	if (colorBy=="--"){
		document.getElementById("coloring").value = "ep";
		colorBy = document.getElementById("coloring").value;
	}
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
		colorScale = cHIColorScale;
		nodes.map(function(d) { d.coloring = d.cHI; });
	}
	else if (colorBy == "date") {
		colorScale = dateColorScale;
		nodes.map(function(d) { d.coloring = d.num_date; });
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
	if (colorBy == "region" || colorBy == "date") {
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
		colorByGenotypePosition(positions_list, 'nuc');
	}
	else {
		d3.select("#coloring").each(colorByTrait);
		gt = parse_gt_string(freqdefault);			
		make_gt_chart(gt);
		document.getElementById("gtspec").value = freqdefault;
	}
}

function colorByGenotypePosition (positions, gene) {
	var gts = nodes.map(function (d) {
		var tmp = [];
		for (var i=0; i<positions.length; i++){
			tmp[tmp.length] = (positions[i]+1)+stateAtPosition(d.clade, gene, positions[i]);
		}
		d.coloring = gene+':'+tmp.join(" / "); 
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

	if ((positions.length==1)&&(filtered_gts.length>1)){
		var tmp_gts=[];
		for (var ii=0; ii<filtered_gts.length; ii+=1){
			tmp_gts.push(["global", filtered_gts[ii]])
		}
		make_gt_chart(tmp_gts);
		document.getElementById("gtspec").value = tmp_gts.map( function (d) {return d[1];}).join(', ');
	}
}

function colorByHIDistance(){
	correctVirus = document.getElementById("virus").checked;
	correctPotency = document.getElementById("serum").checked;
	predictedHI = document.getElementById("HIPrediction").checked;
	colorBy = 'HI_dist'
	if (typeof(focusNode)=="undefined"){
		rootNode.mean_potency=0;
		rootNode.avidity=0;
		focusNode=rootNode;
	}
	treeplot.selectAll(".serum")
		.style("fill", function (d){if (d==focusNode) {return '#FF3300';} else {return '#555555';}})
		.style("font-size", function (d) {if (d==focusNode) {return "32px";} else {return "24px";}})
		.text(function (d) {if (d==focusNode) {return '\uf05b';} else {return '\uf10c';}});

	calcHIpred(focusNode, rootNode);
	calcHImeasured(focusNode, rootNode);

	console.log("Color by HI Distance from "+focusNode.strain);
	console.log("Using predictedHI: "+predictedHI);
	console.log("correcting for virus effect: "+correctVirus);
	console.log("correction for serum effect: "+correctPotency);

	colorScale = HIColorScale;
	if (HIPrediction){
		nodes.map(function(d) { d.coloring = d.HI_dist_pred;});
	}else{
		nodes.map(function(d) { d.coloring = d.HI_dist_meas;});		
	}

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


var genotypeColoringEvent;
d3.select("#gt-color")
	.on("keyup", function(){
		if (typeof genotypeColoringEvent != "undefined"){clearTimeout(genotypeColoringEvent);}
		genotypeColoringEvent = setTimeout(colorByGenotype, 200);
	});



