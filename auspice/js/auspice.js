function gatherTips(node, tips) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			gatherTips(node.children[i], tips);
		}
	}
	else {
		tips.push(node);
	}
	return tips;
}

function gatherInternals(node, internals) {
	if (typeof node.children != "undefined") {
		internals.push(node);
		for (var i=0, c=node.children.length; i<c; i++) {
			gatherInternals(node.children[i], internals);
		}
	}
	return internals;
}

function getVaccines(tips) {
	vaccines = [];
	tips.forEach(function (tip) {
		if (vaccineStrains.indexOf(tip.strain) != -1) {
			tip.choice = vaccineChoice[tip.strain];
			vaccines.push(tip);
		}
	})
	return vaccines;
}

function setDistances(node) {
	if (typeof node.ep == "undefined") {
		node.ep = 0.0;
	}
	if (typeof node.ne == "undefined") {
		node.ne = 0.0;
	}
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			setDistances(node.children[i]);
		}
	}
}

function calcBranchLength(node){
	if (typeof node.children != "undefined") {
	for (var i=0, c=node.children.length; i<c; i++) {
		calcBranchLength(node.children[i]);
		node.children[i].branch_length = node.children[i].xvalue-node.xvalue;
	}
	}
};

/**
sets each node in the tree to alive=true if it has at least one descendent with current=true
**/
function setNodeAlive(node){
	if (typeof node.children != "undefined") {
		var aliveChildren=false;
		for (var i=0, c=node.children.length; i<c; i++) {
			setNodeAlive(node.children[i]);
			aliveChildren = aliveChildren||node.children[i].alive
		}   
		node.alive = aliveChildren;
	}else{
		node.alive = node.current;
	}
};

/**
 * for each node, calculate the exponentially attenuated tree length below the node
 * the polarizer is send "up", i.e. to parents
**/
function calcUpPolarizers(node){
	node.up_polarizer = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
		calcUpPolarizers(node.children[i]);
		node.up_polarizer += node.children[i].up_polarizer;
		}
	}
	bl =  node.branch_length/LBItau;
	node.up_polarizer *= Math.exp(-bl);
	if (node.alive){ // only alive branches contribute anything
		node.up_polarizer += LBItau*(1-Math.exp(-bl));
	}
};

/**
 * for each node, calculate the exponentially attenuated tree length above the node,
 * that is "outside" the clade defined by this node. this down polarizer is send to children
**/
function calcDownPolarizers(node){
	if (typeof node.children != "undefined") {
	for (var i1=0; i1<node.children.length; i1++) {
		node.children[i1].down_polarizer = node.down_polarizer;
		for (var i2=0; i2<node.children.length; i2++) {
			if (i1!=i2){
			node.children[i1].down_polarizer += node.children[i2].up_polarizer;
			}
		}
		// account for the attenuation over the branch_length 
		bl =  node.children[i1].branch_length/LBItau;
		node.children[i1].down_polarizer *= Math.exp(-bl);
		if (node.children[i1].alive) { //the branch contributes only when the node is alive 
			node.children[i1].down_polarizer += LBItau*(1-Math.exp(-bl));
		}
		calcDownPolarizers(node.children[i1]);
	}
	}
};

function calcPolarizers(node){
	calcUpPolarizers(node);
	node.down_polarizer = 0; // set the down polarizer of the root to 0
	calcDownPolarizers(node);
};

/**
 * calculate the LBI for all nodes downstream of node
 * allnodes is provided for easy normalization at the end
**/
function calcLBI(node, allnodes){
	setNodeAlive(node);
	calcPolarizers(node);
	allnodes.forEach(function (d) {
		d.LBI=0;
		d.LBI+=d.down_polarizer;
		if (typeof d.children != "undefined") {
			for (var i=0; i<d.children.length; i++) {
				d.LBI += d.children[i].up_polarizer;
			}
		}
	});
	// normalize the LBI to range [0,1]
	maxLBI = d3.max(allnodes.map(function (d) {return d.LBI;}));
	allnodes.forEach(function (d){ d.LBI /= maxLBI;});
};

/**
 * for each node, calculate the derivative of the frequency tranjectory. if none exists, copy parent
**/
function calcDfreq(node, freq_ii){
	if (typeof node.children != "undefined") {
		for (var i1=0; i1<node.children.length; i1++) {
			if (node.children[i1].freq["global"] != "undefined"){
				var tmp_freq = node.children[i1].freq["global"]
				node.children[i1].dfreq = 0.5*(tmp_freq[freq_ii] - tmp_freq[freq_ii-dfreq_dn])/(tmp_freq[freq_ii] + tmp_freq[freq_ii-dfreq_dn] + 0.1);
			}else{
				node.children[i1].dfreq = node.dfreq;
			}
			calcDfreq(node.children[i1], freq_ii);
		}
	}
};


/**
 * for each node, calculate the number of tips in the currently selected time window. 
**/
function calcTipCounts(node){
	node.tipCount = 0;
	if (typeof node.children != "undefined") {
		for (var i=0; i<node.children.length; i++) {
			calcTipCounts(node.children[i]);
			node.tipCount += node.children[i].tipCount;
		}
	}
	else if (node.current){ 
		node.tipCount = 1;
	}
};

function minimumAttribute(node, attr, min) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			min = minimumAttribute(node.children[i], attr, min);
		}
	}
	else {
		if (node[attr] < min) {
			min = node[attr];
		}
	}
	return min;
}

function maximumAttribute(node, attr, max) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			max = maximumAttribute(node.children[i], attr, max);
		}
	}
	else {
		if (node[attr] > max) {
			max = node[attr];
		}
	}
	return max;
}

function contains(arr, obj) {
    for(var i=0; i<arr.length; i++) {
        if (arr[i] == obj) return true;
    }
}


function branchStrokeColor(col) {
	var modCol = d3.interpolateRgb(col, "#BBB")(0.6);
	return d3.rgb(modCol).toString();
}

function tipStrokeColor(col) {
	return d3.rgb(col).toString();	
}

function tipFillColor(col) {
	return d3.rgb(col).brighter([0.65]).toString();
}

function treePlotHeight(width) {
	return 400 + 0.35*width;
}

var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);

var width = containerWidth,
	height = treePlotHeight(containerWidth);

var cladeToSeq = {}

var globalDate = new Date();
var ymd_format = d3.time.format("%Y-%m-%d");

var LBItau = 0.0008,
	time_window = 1.0;  // layer of one year that is considered current or active

var tree = d3.layout.tree()
	.size([height, width]);

var treeplot = d3.select("#treeplot")
	.attr("width", width)
	.attr("height", height);

var legend = d3.select("#legend")
	.attr("width", 280)
	.attr("height", 100);

var virusTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
	
		string = "";
				
		// safe to assume the following attributes
		if (typeof d.strain != "undefined") {
			string += d.strain;
		}
		string += "<div class=\"smallspacer\"></div>";
		
		string += "<div class=\"smallnote\">";		
		
		// check if vaccine strain
		if (vaccineStrains.indexOf(d.strain) != -1) {
			string += "Vaccine strain<br>";
			var vaccine_date = new Date(vaccineChoice[d.strain]);

			string += "First chosen " + vaccine_date.toLocaleString("en-us", { month: "short" }) + " " + vaccine_date.getFullYear() + "<br>";
			string += "<div class=\"smallspacer\"></div>";
		}			
		
		if (typeof d.country != "undefined") {
			string += d.country.replace(/([A-Z])/g, ' $1');
		}
		if (typeof d.date != "undefined") {
			string += ", " + d.date;
		}
		if ((typeof d.db != "undefined") && (typeof d.accession != "undefined") && (d.db == "GISAID")) {
			string += "<br>GISAID ID: EPI" + d.accession;
		}
		if (typeof d.lab != "undefined") {
			if (d.lab != "") {
				string += "<br>Source: " + d.lab.substring(0,25);
				if (d.lab.length>25) string += '...';
			}
		}			
		string += "</div>";
		
		string += "<div class=\"smallspacer\"></div>";
				
		// following may or may not be present
		string += "<div class=\"smallnote\">";
		if (typeof d.ep != "undefined") {
			string += "Epitope distance: " + d.ep + "<br>";
		}
		if (typeof d.ne != "undefined") {
			string += "Non-epitope distance: " + d.ne + "<br>";
		}
		if (typeof d.rb != "undefined") {
			string += "Receptor binding distance: " + d.rb + "<br>";
		}
		if (typeof d.LBI != "undefined") {
			string += "Local branching index: " + d.LBI.toFixed(3) + "<br>";
		}
		string += "</div>";
		return string;
	});
treeplot.call(virusTooltip);

var linkTooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 12])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			string += "Frequency: " + (100 * d.frequency).toFixed(1) + "%"
			if (d.aa_muts.length){
				string+="<br>Mutations: "+d.aa_muts.replace(/,/g, ', ');
			}
		}
		return string;
	});
treeplot.call(linkTooltip);

width = parseInt(d3.select(".freqplot-container").style("width"), 10);
var position = "right";
if (width < 600) {
	position = "bottom";
}

var gt_chart = c3.generate({
	bindto: '#gtchart',
	size: {width: width, height: 350},
	legend: {position: position},
  	color: {
        pattern: ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]
    },
	axis: {
		y: {
			label: {
				text: 'frequency',
				position: 'outer-middle'	
			},
			tick: {
				values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
				outer: false
			},
            min: 0,			
			max: 1	
		},
		x: {
			label: {
				text: 'time',
				position: 'outer-center'	
			},
			tick: {
				values: time_ticks,
				outer: false				
			}
		}
	},			
	data: {
		x: 'x',
		columns: []
	}
});


d3.json("/data/" + file_prefix + "tree.json", function(error, root) {

	if (error) return console.warn(error);

	var nodes = tree.nodes(root),
		links = tree.links(nodes);
	var tree_legend;
	var rootNode = nodes[0];
	if (typeof rootNode.pivots != "undefined"){
		var dt = rootNode.pivots[1]-rootNode.pivots[0];		
	}else{
		var dt = 1.0/12;
	}
	var tips = gatherTips(rootNode, []);
	var internals = gatherInternals(rootNode, []);
	calcBranchLength(rootNode);
	rootNode.branch_length= 0.01;
	nodes.forEach(function (d) {d.dateval = new Date(d.date)});

	var vaccines = getVaccines(tips);

	var xValues = nodes.map(function(d) {
		return +d.xvalue;
	});

	var yValues = nodes.map(function(d) {
		return +d.yvalue;
	});

	var dateValues = nodes.filter(function(d) {
		return typeof d.date === 'string';
		}).map(function(d) {
		return new Date(d.date);
	});

	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)])
		.range([10, width-10]);

	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)])
		.range([10, height-10]);

	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);
	});

	var earliestDate = new Date(d3.min(dateValues));
	earliestDate.setDate(earliestDate.getDate() + 180);

	var dateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true]);	

	var niceDateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true])
		.nice(d3.time.month);

	var recencySizeScale = d3.scale.threshold()
		.domain([0.0, 1.0])
		.range([0, 4, 0]);

	var recencyVaccineSizeScale = d3.scale.threshold()
		.domain([0.0])
		.range([0, 8]);

	var recencyLinksSizeScale = d3.scale.threshold()
		.domain([0.0])
		.range([0, 2]);

	var colors = ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"];
	var colorBy = document.getElementById("coloring").value;
	
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
		.domain(([-1.0, -0.8, -0.6,-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]).map(function(d){return Math.round(d*dt*dfreq_dn*100)/100;}))
		.range(colors);

	var colorScale;
	
	var freqScale = d3.scale.linear()
		.domain([0, 1])
		.range([1.5, 4.5]);

	var regions = ["Africa", "SouthAmerica", "WestAsia", "Oceania", "Europe", "JapanKorea", "NorthAmerica", "SoutheastAsia", "India", "China"]
	var regionColors = ["#5097BA", "#60AA9E", "#75B681", "#8EBC66", "#AABD52", "#C4B945", "#D9AD3D", "#E59637", "#E67030", "#DF4327"]

	var regionColorScale = d3.scale.ordinal()
		.domain(regions)
		.range(regionColors);

	var genotypeColors = ["#60AA9E", "#D9AD3D", "#5097BA", "#E67030", "#8EBC66", "#E59637", "#AABD52", "#DF4327", "#C4B945", "#75B681"]

	function calcNodeAges(tw){
		tips.forEach(function (d) {
			var date = new Date(d.date);
			var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
			var diffYears = (globalDate.getTime() - date.getTime()) / oneYear;
			d.diff = diffYears;
			if (d.diff > 0 && d.diff < tw){
				d.current  = true;
			}else{
				d.current = false;
			}
		});
	};
	
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
	
	function adjust_freq_by_date() {
		calcTipCounts(rootNode);
		var tipCount = rootNode.tipCount;		
		console.log("Total tipcount: " + tipCount);	
		nodes.forEach(function (d) {
			d.frequency = (d.tipCount)/tipCount;
		});
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
			.style("stroke", function(d) {
					if (colorScale != regionColorScale) {
						var col = colorScale(d.target.adj_coloring);
					}
					else {
						var col = "#AAA";
					}
					return branchStrokeColor(col);
				});
			
		d3.selectAll(".tip")
			.attr("r", function(d) {
				return recencySizeScale(d.diff);
			})
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

    var legendRectSize = 15;
    var legendSpacing = 4;
    function makeLegend(){
    	
    	d3.select("#legend-title").text(function(d){
    		if (colorBy == "ep") {
    			return "Epitope mutations"
    		}
    		if (colorBy == "ne") {
    			return "Non-epitope mutations"
    		}
    		if (colorBy == "rb") {
    			return "Receptor binding mutations"
    		}
    		if (colorBy == "lbi") {
    			return "Local branching index"
    		}
   			if (colorBy == "region") {
    			return "Region"
    		}
   			if (colorBy == "genotype") {
    			return "Genotype"
    		}
   			if (colorBy == "dfreq") {
   				var tmp_nmonth = Math.round(12*dfreq_dn*dt);
   				var tmp_text = "Freq. change ("+tmp_nmonth+" month";
   				if (tmp_nmonth>1){
   					tmp_text+='s';
   				}
    			return tmp_text+')';
    		}
    	});
    
		var tmp_leg = legend.selectAll(".legend")
		    .data(colorScale.domain())
		    .enter().append('g')
		    .attr('class', 'legend')
		    .attr('transform', function(d, i) {
		    	var stack = 5;
				var height = legendRectSize + legendSpacing;
				var fromRight = Math.floor(i / stack);
				var fromTop = i % stack;
				var horz = fromRight * 145 + 5;				
				var vert = fromTop * height + 5;
				return 'translate(' + horz + ',' + vert + ')';
		    	});
		tmp_leg.append('rect')
		    .attr('width', legendRectSize)
		    .attr('height', legendRectSize)
		    .style('fill', function (d) {
		   		var col = colorScale(d);
		   		return d3.rgb(col).brighter([0.35]).toString();
		    })
		    .style('stroke', function (d) {
		    	var col = colorScale(d);
		    	return tipStrokeColor(col);
		    });
		
		tmp_leg.append('text')
		    .attr('x', legendRectSize + legendSpacing + 5)
		    .attr('y', legendRectSize - legendSpacing)
		    .text(function(d) {
		    	return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
		    });		
		return tmp_leg;
    }

    function removeLegend(){
    	legend.selectAll('.legend')
    		.remove();
    }

	calcNodeAges(time_window);
	calcLBI(rootNode, nodes, false);
	calcDfreq(rootNode, freq_ii);
	var freq_ii = rootNode.pivots.length - 1;
	console.log(freq_ii);
	colorByTrait();
	adjust_coloring_by_date();
	adjust_freq_by_date();

	var link = treeplot.selectAll(".link")
		.data(links)
		.enter().append("polyline")
		.attr("class", "link")
		.attr("points", function(d) {
			var mod = 0.5 * freqScale(d.target.frequency) - freqScale(0);
			return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
			+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
			+ (d.target.x).toString() + "," + d.target.y.toString()
		})
		.style("stroke-width", function(d) {
			return freqScale(d.target.frequency);
		})
		.style("stroke", function(d) {
				var col = colorScale(d.target.adj_coloring);
				return branchStrokeColor(col);
			})		
		.style("cursor", "pointer")
		.on('mouseover', function(d) {
			linkTooltip.show(d.target, this);
			var plot_data = [['x'].concat(rootNode["pivots"])];
			var reg = "global";
			if (d.target.freq[reg] != "undefined"){
				plot_data[plot_data.length] = [reg].concat(d.target.freq[reg]);				
			}
			if (plot_data.length > 1) {
				if (plot_data[1][0] == "global") {
					plot_data[1][0] = "clade";
				}
			}
			gt_chart.load({
		       	columns: plot_data
			});
		})
		.on('mouseout', linkTooltip.hide)		
		.on('click', function(d) {
			var dMin = minimumAttribute(d.target, "xvalue", d.target.xvalue),
				dMax = maximumAttribute(d.target, "xvalue", d.target.xvalue),
				lMin = minimumAttribute(d.target, "yvalue", d.target.yvalue),
				lMax = maximumAttribute(d.target, "yvalue", d.target.yvalue);
			if (dMax > dMin && lMax > lMin) {
				rescale(dMin, dMax, lMin, lMax);
			}
			else {
				dMin = minimumAttribute(d.source, "xvalue", d.source.xvalue),
				dMax = maximumAttribute(d.source, "xvalue", d.source.xvalue),
				lMin = minimumAttribute(d.source, "yvalue", d.source.yvalue),
				lMax = maximumAttribute(d.source, "yvalue", d.source.yvalue);			
				rescale(dMin, dMax, lMin, lMax);
			}
		});

	var tipCircles = treeplot.selectAll(".tip")
		.data(tips)
		.enter()
		.append("circle")
		.attr("class", "tip")
		.attr("id", function(d) { return (d.strain).replace(/\//g, ""); })
		.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; })
		.attr("r", function(d) {
			return recencySizeScale(d.diff);
		})
		.style("fill", function(d) {
			var col = colorScale(d.adj_coloring);
			return tipFillColor(col);
		})
		.style("stroke", function(d) {
			var col = colorScale(d.adj_coloring);
			return tipStrokeColor(col);
		})
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('click', function(d) {
			if ((typeof d.db != "undefined") && (d.db == "GISAID") && (typeof d.accession != "undefined")) {
				var url = "http://gisaid.org/EPI/"+d.accession;
				console.log("opening url "+url);
				var win = window.open(url, '_blank');
  				win.focus();
  			}	
  		})		
		.on('mouseout', virusTooltip.hide);

	var vaccineCircles = treeplot.selectAll(".vaccine")
		.data(vaccines)
		.enter()
		.append("text")
		.attr("class", "vaccine")
		.attr("x", function(d) {return d.x})
		.attr("y", function(d) {return d.y})
		.attr('text-anchor', 'middle')
		.attr('dominant-baseline', 'central')
		.style("font-size", "28px")
		.style('font-family', 'FontAwesome')
		.style("fill", "#555555")
		.text(function(d) { return '\uf00d'; })
		.style("cursor", "default")
		.on('mouseover', function(d) {
			virusTooltip.show(d, this);
		})
		.on('mouseout', virusTooltip.hide);


	var drag = d3.behavior.drag()
		.origin(function(d) { return d; })
		.on("drag", dragged)
		.on("dragstart", function() {
			d3.selectAll(".date-input-text").style("fill", "#5DA8A3");
			d3.selectAll(".date-input-marker").style("fill", "#5DA8A3");
		})
		.on("dragend", function() {
			d3.selectAll(".date-input-text").style("fill", "#CCC");
			d3.selectAll(".date-input-marker").style("fill", "#CCC");
			dragend();
		});

	function dragged(d) {

		d.date = dateScale.invert(d3.event.x);
		d.x = dateScale(d.date);
		d3.selectAll(".date-input-text")
			.attr("dx", function(d) {return 0.5*d.x})
			.text(function(d) {
				var format = d3.time.format("%Y %b %-d");
				return format(d.date)
			});
		d3.selectAll(".date-input-marker")
			.attr("cx", function(d) {return d.x});
		globalDate = d.date;

		calcNodeAges(time_window);
		treeplot.selectAll(".link")
			.style("stroke", function(d){return "#ccc";})

		treeplot.selectAll(".tip")
			.attr("r", function(d) {
				return recencySizeScale(d.diff);
			})
			.style("fill", "#CCC")
			.style("stroke", "#AAA");

		treeplot.selectAll(".vaccine")
			.style("visibility", function(d) {
				var date = new Date(d.choice);
				var oneYear = 365.25*24*60*60*1000; // days*hours*minutes*seconds*milliseconds
				var diffYears = (globalDate.getTime() - date.getTime()) / oneYear;
				if (diffYears > 0) { return "visible"; }
				else { return "hidden"; }
			});					
		
	}

	function dragend() {
		var num_date = globalDate/1000/3600/24/365.25+1970;	
		for (var ii=0; ii<rootNode.pivots.length-1; ii++){
			if (rootNode.pivots[ii]<num_date && rootNode.pivots[ii+1]>=num_date){
				freq_ii=Math.max(dfreq_dn,ii+1);
			}
		}
		console.log("changed frequency index to "+freq_ii+" date cut off is "+num_date);
		console.log("recalculating node ages");
		calcNodeAges(time_window);
		console.log("adjusting node colors");
		adjust_coloring_by_date();
		console.log("updating frequencies");
		adjust_freq_by_date();
		
		if (colorBy == "genotype") {
			colorByGenotype();
		}

		if (colorBy!="genotype"){
			d3.selectAll(".link")
				.transition().duration(500)
				.attr("points", function(d) {
					var mod = 0.5 * freqScale(d.target.frequency) - freqScale(0);				
					return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
					+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
					+ (d.target.x).toString() + "," + d.target.y.toString()
				})
				.style("stroke-width", function(d) {
					return freqScale(d.target.frequency);
				})
				.style("stroke", function(d) {
					if (colorScale != regionColorScale) {
						var col = colorScale(d.target.adj_coloring);
					}
					else {
						var col = "#AAA";
					}
					return branchStrokeColor(col);
				});				
				
			d3.selectAll(".tip")
				.transition().duration(500)
				.attr("r", function(d) {
					return recencySizeScale(d.diff);
				})
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
		}
	}

	var counterData = {}
	counterData['date'] = globalDate
	counterData['x'] = dateScale(globalDate)

	d3.select("#date-input")
		.attr("width", 240)
		.attr("height", 65);

	var counter = d3.select("#date-input").selectAll(".date-input-text")
		.data([counterData])
		.enter()
		.append("text")
		.attr("class", "date-input-text")
		.attr("text-anchor", "left")
		.attr("dx", function(d) {return 0.5*d.x})		
		.attr("dy", "1.0em")
		.text(function(d) {
			var format = d3.time.format("%Y %b %-d");
			return format(d.date)
		})
		.style("cursor", "pointer")
		.call(drag);

	var customTimeFormat = d3.time.format.multi([
		[".%L", function(d) { return d.getMilliseconds(); }],
		[":%S", function(d) { return d.getSeconds(); }],
		["%I:%M", function(d) { return d.getMinutes(); }],
		["%I %p", function(d) { return d.getHours(); }],
		["%a %d", function(d) { return d.getDay() && d.getDate() != 1; }],
		["%b %d", function(d) { return d.getDate() != 1; }],
		["%b", function(d) { return d.getMonth(); }],
		["%Y", function() { return true; }]
	]);

	var dateAxis = d3.svg.axis()
		.scale(niceDateScale)
		.orient('bottom')
		.ticks(5)
		.tickFormat(customTimeFormat)
		.outerTickSize(2)
		.tickPadding(8);

	d3.select("#date-input").selectAll(".date-input-axis")
		.data([counterData])
		.enter()
		.append("g")
		.attr("class", "date-input-axis")
		.attr("transform", "translate(0,35)")
		.call(dateAxis);

	var marker = d3.select("#date-input").selectAll(".date-input-marker")
		.data([counterData])
		.enter()
		.append("circle")
		.attr("class", "date-input-marker")
		.attr("cx", function(d) {return d.x})
		.attr("cy", 35)
		.attr("r", 5)
		.style("fill", "#CCC")
		.style("stroke", "#777")		
		.style("cursor", "pointer")
		.call(drag);

	d3.select("#reset")
		.on("click", function(d) {
			var dMin = d3.min(xValues),
				dMax = d3.max(xValues),
				lMin = d3.min(yValues),
				lMax = d3.max(yValues);
			rescale(dMin, dMax, lMin, lMax);
		})

	function rescale(dMin, dMax, lMin, lMax) {

		var speed = 1500;
		xScale.domain([dMin,dMax]);
		yScale.domain([lMin,lMax]);

		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});

		treeplot.selectAll(".tip").data(tips)
			.transition().duration(speed)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".vaccine").data(vaccines)
			.transition().duration(speed)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });

		treeplot.selectAll(".internal").data(internals)
			.transition().duration(speed)
			.attr("x", function(d) {
				if (typeof d.frequency != "undefined") {
					return d.x - 5*Math.sqrt(d.frequency) - 0.5;
				}
				else {
					return d.x - 1;
				}
			})
			.attr("y", function(d) {
				if (typeof d.frequency != "undefined") {
					return d.y - 5*Math.sqrt(d.frequency) - 0.5;
				}
				else {
					return d.y - 1;
				}
			});

		treeplot.selectAll(".link").data(links)
			.transition().duration(speed)
			.attr("points", function(d) {
				return (d.source.x).toString() + "," + d.source.y.toString() + " "
				+ (d.source.x).toString() + "," + d.target.y.toString() + " "
				+ (d.target.x).toString() + "," + d.target.y.toString()
			});
			
		treeplot.selectAll(".annotation").data(clades)
			.transition().duration(speed)
			.attr("x", function(d) {
				return xScale(d[1]) - 8;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 8;
			});			

	}	

	d3.select(window).on('resize', resize); 
	
	function resize() {
	
		var containerWidth = parseInt(d3.select(".treeplot-container").style("width"), 10);
		var width = containerWidth,
			height = treePlotHeight(containerWidth);
			
		d3.select("#treeplot")
			.attr("width", width)
			.attr("height", height);			
			
		xScale.range([10, width-10]);
		yScale.range([10, height-10]);
		
		nodes.forEach(function (d) {
			d.x = xScale(d.xvalue);
			d.y = yScale(d.yvalue);
		});		
		
		treeplot.selectAll(".tip").data(tips)
			.attr("cx", function(d) { return d.x; })
			.attr("cy", function(d) { return d.y; });

		treeplot.selectAll(".vaccine").data(vaccines)
			.attr("x", function(d) { return d.x; })
			.attr("y", function(d) { return d.y; });
			
		treeplot.selectAll(".link").data(links)
			.attr("points", function(d) {
				return (d.source.x).toString() + "," + d.source.y.toString() + " "
				+ (d.source.x).toString() + "," + d.target.y.toString() + " "
				+ (d.target.x).toString() + "," + d.target.y.toString()
			});
			
		treeplot.selectAll(".annotation").data(clades)
			.attr("x", function(d) {
				return xScale(d[1]) - 8;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 8;
			});

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

	function onSelect(tip) {
		d3.select("#"+(tip.strain).replace(/\//g, ""))
			.call(function(d) {
				virusTooltip.show(tip, d[0][0]);
			});
	}

	var mc = autocomplete(document.getElementById('search'))
		.keys(tips)
		.dataField("strain")
		.placeHolder("search strains...")
		.width(800)
		.height(500)
		.onSelected(onSelect)
		.render();



	d3.select("#gt-color")
		.on("keyup", colorByGenotype);

	tree_legend = makeLegend();

	// add clade labels
	clades = rootNode["clade_annotations"];
	console.log(clades);
	var clade_annotations = treeplot.selectAll('.annotation')
		.data(clades)
		.enter()
		.append("text")
		.attr("class", "annotation")
		.attr("x", function(d) {
			return xScale(d[1]) - 8;
		})
		.attr("y", function(d) {
			return yScale(d[2]) - 8;
		})
		.style("text-anchor", "end")
		.text(function (d) {
			return d[0];
		});


});

d3.json("/data/" + file_prefix + "meta.json", function(error, json) {
	if (error) return console.warn(error);
	d3.select("#updated").text(json['updated']);
	commit_id = json['commit'];
	short_id = commit_id.substring(0, 6);	
	d3.select("#commit")
		.append("a")
		.attr("href", "http://github.com/blab/nextflu/commit/" + commit_id)
		.text(short_id);

});

d3.json("/data/" + file_prefix + "sequences.json", function(error, json) {
	if (error) return console.warn(error);
	cladeToSeq=json;
});

d3.json("/data/" + file_prefix + "frequencies.json", function(error, json){
	console.log(error);
	var pivots= json["mutations"]["global"]["pivots"].map(function (d) {return Math.round(parseFloat(d)*100)/100;});
	var ticks = [Math.round(pivots[0])];
	var step = Math.round((pivots[pivots.length-1]-pivots[0])/6*10)/10;
	while (ticks[ticks.length-1]<pivots[pivots.length-1]){
		ticks.push(Math.round((ticks[ticks.length-1]+step)*10)/10);
	}
	//gt_chart.axis.x.values = ticks;
	/**
		parses a genotype string into region and positions
	**/
	function parse_gt_string(gt){
		separate_plots = gt.split(',');
		mutations = separate_plots.map(
			function (d) {	var tmp = d.split(/[\s//]/); //FIXME: make more inclusive
							var region;
							var positions = [];
							for (var i=0; i<tmp.length; i++){
								if (contains(["EU","NA","AS","OC"], tmp[i])){
									region = tmp[i];
								}else{
									if (tmp[i].length>0) positions.push(tmp[i]);
								}
							}
							if (typeof region == "undefined") region="global";
							// sort of this is a multi mutation genotype
							if (positions.length>1){
								positions.sort(function (a,b){
									return parseInt(a.substring(0,a.length-1)) - parseInt(b.substring(0,b.length-1));
								});
							}
							return [region, positions.join('/')];});
		return mutations;
	};

	/**
	loops over all genotypes from a certain region and sums the frequency contributions
	of the genotype matches at the specified positions
	**/
	function get_frequencies(region, gt){
		var freq = [];
		for (var pi=0; pi<pivots.length; pi++){freq[freq.length]=0;}
		if (json["clades"][region][gt.toLowerCase()]!=undefined) {
			console.log(gt+" found as clade");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["clades"][region][gt.toLowerCase()][pi];
			}
		}
		else if ((typeof json["genotypes"] !="undefined") && (json["genotypes"][region][gt]!=undefined)) {
			console.log(gt+" found as genotype");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["genotypes"][region][gt][pi];
			}
		}else if (json["mutations"][region][gt]!=undefined) {
			console.log(gt+" found as mutation");
			for (var pi=0; pi<freq.length; pi++){
				freq[pi]+=json["mutations"][region][gt][pi];
			}
		}
		return freq.map(function (d) {return Math.round(d*100)/100;});
	};

	function make_gt_chart(gt){
		var tmp_data = [];
		var tmp_trace = ['x'];
		tmp_data.push(tmp_trace.concat(pivots));
		gt.forEach(function (d) {
			var region = d[0];
			var genotype = d[1];
			var freq = get_frequencies(region, genotype);
			var tmp_trace = genotype.toString().replace(/,/g, ', ');
			if (region != "global") {
				tmp_trace = region + ':\t' + tmp_trace;
			}
			tmp_data.push([tmp_trace].concat(freq));
		});
		gt_chart.load({
	       	columns: tmp_data,
	       	unload: true
		});
	}

	d3.select("#plotfreq")
		.on("click", function (){
			gt = parse_gt_string(document.getElementById("gtspec").value);			
			make_gt_chart(gt);
		});
	make_gt_chart(parse_gt_string(document.getElementById("gtspec").value));
});
