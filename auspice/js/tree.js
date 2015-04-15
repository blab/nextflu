var recencySizeScale = d3.scale.threshold()
	.domain([0.0, time_window])
	.range([0, 4, 0]);

var recencyVaccineSizeScale = d3.scale.threshold()
	.domain([0.0])
	.range([0, 8]);

var recencyLinksSizeScale = d3.scale.threshold()
	.domain([0.0])
	.range([0, 2]);

var freqScale = d3.scale.linear()
	.domain([0, 1])
	.range([1.5, 4.5]);

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

function branchPoints(d) {
	var mod = 0.5 * freqScale(d.target.frequency) - freqScale(0);
	return (d.source.x-mod).toString() + "," + d.source.y.toString() + " "
		+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
		+ (d.target.x).toString() + "," + d.target.y.toString();
}

function branchStrokeWidth(d) {
	return freqScale(d.target.frequency);
}

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
		.attr("points", branchPoints)
		.style("stroke-width", branchStrokeWidth)
		.style("stroke", branchStrokeColor)		
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
		.attr("r", tipRadius)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor)
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
			.attr("points", branchPoints);
			
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
			.attr("points", branchPoints);
			
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
			.attr("r", tipRadius);			
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

