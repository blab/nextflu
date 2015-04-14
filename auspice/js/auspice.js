

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


var recencySizeScale = d3.scale.threshold()
	.domain([0.0, time_window])
	.range([0, 4, 0]);


var recencyVaccineSizeScale = d3.scale.threshold()
	.domain([0.0])
	.range([0, 8]);

var recencyLinksSizeScale = d3.scale.threshold()
	.domain([0.0])
	.range([0, 2]);

function tipRadius(d) {
	var radius = 0;
	if (d.region == restrictTo || restrictTo == "all") {
		radius = recencySizeScale(d.diff);
	}
	else {
		radius = 0;
	}
	return radius;
}


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

function tree_init(){
	calcBranchLength(rootNode);
	rootNode.branch_length= 0.01;	
	rootNode.dfreq = 0.0;
	freq_ii = 1;
	if (typeof rootNode.pivots != "undefined"){
		dt = rootNode.pivots[1]-rootNode.pivots[0];		
	}else{
		dt = 1.0/12;
	}
	if (typeof rootNode.pivots != "undefined") {
		if (typeof rootNode.pivots.length != "undefined") {
			freq_ii = rootNode.pivots.length - 1;
		}
	}
	dfreqColorScale = d3.scale.linear()
		.domain(([-1.0, -0.8, -0.6,-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8]).map(function(d){return Math.round(d*dt*dfreq_dn*100)/100;}))
		.range(colors)

	calcNodeAges(LBItime_window);
	calcLBI(rootNode, nodes, false);
	calcDfreq(rootNode, freq_ii);
	colorByTrait();
	adjust_coloring_by_date();
	adjust_freq_by_date();
	tree_legend = makeLegend();
}

d3.json("/data/" + file_prefix + "tree.json", function(error, root) {

	if (error) return console.warn(error);

	nodes = tree.nodes(root);
	links = tree.links(nodes);
	var tree_legend;
	rootNode = nodes[0];
	tips = gatherTips(rootNode, []);
	internals = gatherInternals(rootNode, []);
	vaccines = getVaccines(tips);

	var xValues = nodes.map(function(d) {
		return +d.xvalue;
	});

	var yValues = nodes.map(function(d) {
		return +d.yvalue;
	});

	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)])
		.range([10, treeWidth-10]);

	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)])
		.range([10, treeHeight-10]);

	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);
	});

	date_init();
	tree_init();

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
			var dMin = 0.5 * (minimumAttribute(d.target, "xvalue", d.target.xvalue) + minimumAttribute(d.source, "xvalue", d.source.xvalue)),
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
		.attr("r", function(d) { return tipRadius(d); })
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
				return xScale(d[1]) - 6;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 6;
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
				return xScale(d[1]) - 6;
			})
			.attr("y", function(d) {
				return yScale(d[2]) - 6;
			});

	}
	
		
	restrictTo = document.getElementById("region").value;

	function restrictToRegion() {
		restrictTo = document.getElementById("region").value;
		console.log(restrictTo);	
		d3.selectAll(".tip")
			.attr("r", function(d) { return tipRadius(d); });			
	}

	d3.select("#region")
		.style("cursor", "pointer")
		.on("change", restrictToRegion);		

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


	// add clade labels
	clades = rootNode["clade_annotations"];
	console.log(clades);
	var clade_annotations = treeplot.selectAll('.annotation')
		.data(clades)
		.enter()
		.append("text")
		.attr("class", "annotation")
		.attr("x", function(d) {
			return xScale(d[1]) - 6;
		})
		.attr("y", function(d) {
			return yScale(d[2]) - 6;
		})
		.style("text-anchor", "end")
		.text(function (d) {
			return d[0];
		});


});

d3.json("/data/" + file_prefix + "sequences.json", function(error, json) {
	if (error) return console.warn(error);
	cladeToSeq=json;
});

