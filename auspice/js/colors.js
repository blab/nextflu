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

var dfreqColorScale; //defined after loading the tree.


var regionColorScale = d3.scale.ordinal()
	.domain(regions)
	.range(regionColors);

function tipStrokeColor(col) {
	return d3.rgb(col).toString();	
}

function tipFillColor(col) {
	return d3.rgb(col).brighter([0.65]).toString();
}

function getMeanColoring() {	
	var mean = 0;
	var recent_tip_count = 0;
	tips.forEach(function (d) {
		if (d.current) {
			mean += d.coloring;
			recent_tip_count += 1;
		}
	});
	mean = mean / recent_tip_count;
	return mean;
}	

function adjust_coloring_by_date() {
	if (colorBy == "ep" || colorBy == "ne" || colorBy == "rb") {
		var mean = getMeanColoring();
		nodes.forEach(function (d) {
			d.adj_coloring = d.coloring; // - mean;
		});
	}
	else if (colorBy == "lbi") {
		calcLBI(rootNode, nodes, false);
		nodes.forEach(function (d) {
			d.adj_coloring = d.LBI;
		});
	}
	else if (colorBy == "dfreq") {
		calcDfreq(rootNode, freq_ii);
		nodes.forEach(function (d) {
			d.adj_coloring = d.dfreq;
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
		nodes.map(function(d) { d.adj_coloring = d.LBI; });
	}
	else if (colorBy == "dfreq") {
		colorScale = dfreqColorScale;
		nodes.map(function(d) { d.adj_coloring = d.dfreq; });
	}
	else if (colorBy == "region") {
		colorScale = regionColorScale;
	}

	adjust_coloring_by_date();

	treeplot.selectAll(".link")
		.style("stroke", branchStrokeColor);
		
	d3.selectAll(".tip")
		.attr("r", function(d) { return tipRadius(d); })
		.style("fill", function(d) {
			if (colorScale != regionColorScale) {
				var col = colorScale(d.adj_coloring);
			}
			else {
				var col = colorScale(d.region);
			}
			return tipFillColor(col);
		})
		.style("stroke", function(d) {
			if (colorScale != regionColorScale) {
				var col = colorScale(d.adj_coloring);
			}
			else {
				var col = colorScale(d.region);
			}
			return tipStrokeColor(col);
		});
		
	if (typeof tree_legend != undefined){
		removeLegend();
	}
	tree_legend = makeLegend();	 				
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
	var gts = nodes.map(function (d) {var tmp = [];
										for (var i=0; i<positions.length; i++){
											var aa = cladeToSeq[d.clade];
											tmp[tmp.length] = (positions[i]+1)+aa[positions[i]];
										}
										d.color_gt = tmp.join(" / "); 
										return d.color_gt;});
	var unique_gts = d3.set(gts).values();
	var gt_counts = {};
	for (var i=0; i<unique_gts.length; i++){gt_counts[unique_gts[i]]=0;}
	gts.forEach(function (d) {gt_counts[d]+=1;});
	var filtered_gts = unique_gts.filter(function (d) {return gt_counts[d]>=10;});
	filtered_gts.sort(function (a,b){var res;
		if (gt_counts[a]>gt_counts[b]){ res=-1;}
		else if (gt_counts[a]<gt_counts[b]){ res=1;}
		else {res=0;}
		return res;});
	console.log("genotypes passed filtering:"+filtered_gts);
	colorScale = d3.scale.ordinal()
		.domain(filtered_gts)
		.range(genotypeColors);			
		treeplot.selectAll(".link")
		.style("stroke", function(d) {
			var col = colorScale(d.target.color_gt);
			return branchStrokeColor(col);
		});
		treeplot.selectAll(".tip")
		.style("fill", function(d) {
			var col = colorScale(d.color_gt);
			return tipFillColor(col);
		})
		.style("stroke", function(d) {
			var col = colorScale(d.color_gt);
			return tipStrokeColor(col);
		});
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



