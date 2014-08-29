function depthFirstSearch(node) {
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			depthFirstSearch(node.children[i]);
		}
	}
}

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

function setFrequencies(node) {
	if (typeof node.frequency == "undefined") {
		node.frequency = 0.01;
	}
	if (typeof node.children != "undefined") {
		for (var i=0, c=node.children.length; i<c; i++) {
			setFrequencies(node.children[i]);
		}
	}
}

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

var width = 800,
	height = 600;

var tree = d3.layout.tree()
	.size([height, width]);

var svg = d3.select("svg")
	.attr("width", width)
	.attr("height", height);
		
var tooltip = d3.tip()
	.direction('e')
	.attr('class', 'd3-tip')
	.offset([0, 10])
	.html(function(d) {
		string = ""
		if (typeof d.frequency != "undefined") {
			if (d.frequency > 0.01) {
				string = d.frequency;
			}
		}
		if (typeof d.target != "undefined") {
			if (typeof d.target.frequency != "undefined") {
				if (d.target.frequency > 0.01) {
					string = d.target.frequency;
				}
			}	
		}	
		if (typeof d.strain != "undefined") {
			string = d.strain;
		}		
		return string;
	});
	
svg.call(tooltip);		

function rescale(dMin, dMax, lMin, lMax, xScale, yScale, nodes, links, tips, internals) {

	var speed = 1500;
	xScale.domain([dMin,dMax]);
	yScale.domain([lMin,lMax]);
	
	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);			 
	});	
		
	svg.selectAll(".tip").data(tips)
    	.transition().duration(speed)
    	.attr("cx", function(d) { return d.x; })
    	.attr("cy", function(d) { return d.y; }); 
    	
	svg.selectAll(".internal").data(internals)
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
	    
	svg.selectAll(".link").data(links)
    	.transition().duration(speed)
		.attr("points", function(d) {
			var mod = 5*Math.sqrt(d.target.frequency)+0.5;
			return (d.source.x-mod).toString() + "," + d.source.y.toString() + " " 
			+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
			+ (d.target.x).toString() + "," + d.target.y.toString()
		});	   		
		
}

d3.json("https://s3.amazonaws.com/trvrb-augur/tree.json", function(error, root) {
	var nodes = tree.nodes(root),
		links = tree.links(nodes);
	
	var rootNode = nodes[0];
	var tips = gatherTips(rootNode, []);
	var internals = gatherInternals(rootNode, []);
	setFrequencies(rootNode);
		
	var	xValues = nodes.map(function(d) {
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
  	
  	var yearValues = nodes.filter(function(d) {
		return typeof d.date === 'string';
  		}).map(function(d) {
  		return (new Date(d.date)).getFullYear();
  	}); 	  	

	var xScale = d3.scale.linear()
		.domain([d3.min(xValues), d3.max(xValues)])
		.range([10, width-10]);
		
	var yScale = d3.scale.linear()
		.domain([d3.min(yValues), d3.max(yValues)])
		.range([10, height-10]);	
						
	var dateScale = d3.time.scale()
		.domain([d3.min(dateValues), d3.max(dateValues)])
		.range([0, 1]);
		
	var yearScale = d3.scale.ordinal()
		.domain([2014, "undefined", 2011, 2012, 2013])
		.range(["#ff7f0e", "#1f77b4", "#7f7f7f", "#7f7f7f", "#7f7f7f"]);
						
	nodes.forEach(function (d) {
		d.x = xScale(d.xvalue);
		d.y = yScale(d.yvalue);			 
	});

	// straight links
/*	var link = svg.selectAll(".link")
		.data(links)
		.enter().append("line")
		.attr("class", "link")
		.attr("x1", function(d) { return d.source.y; })
	    .attr("y1", function(d) { return d.source.x; })
	    .attr("x2", function(d) { return d.target.y; })
	    .attr("y2", function(d) { return d.target.x; }); */
	    	    	    
	var link = svg.selectAll(".link")
		.data(links)
		.enter().append("polyline")
		.attr("class", "link")
		.attr("points", function(d) {
			var mod = 5*Math.sqrt(d.target.frequency)+0.5;
			return (d.source.x-mod).toString() + "," + d.source.y.toString() + " " 
			+ (d.source.x-mod).toString() + "," + d.target.y.toString() + " "
			+ (d.target.x).toString() + "," + d.target.y.toString()
		})
	    .style("stroke-width", function(d) {
			return 10*Math.sqrt(d.target.frequency)+1;		    	
	    })
		.on('mouseover', tooltip.show)
      	.on('mouseout', tooltip.hide)	    
     	.on('click', function(d) { 
      		var dMin = minimumAttribute(d.target, "xvalue", d.target.xvalue),
      			dMax = maximumAttribute(d.target, "xvalue", d.target.xvalue),
      			lMin = minimumAttribute(d.target, "yvalue", d.target.yvalue),
      			lMax = maximumAttribute(d.target, "yvalue", d.target.yvalue);
      		rescale(dMin, dMax, lMin, lMax, xScale, yScale, nodes, links, tips, internals)
      	});   	    
	    
	var tipCircles = svg.selectAll(".tip")
		.data(tips)
		.enter()
		.append("circle")
		.attr("class", "tip")
		.attr("cx", function(d) {return d.x})
		.attr("cy", function(d) {return d.y})
		.attr("r", 1.5)		
		.style("fill", function(d) { 
			var today = new Date();
			var date = new Date(d.date);
			var diff = (today - date) / (1000*60*60*24*365.25) // convert from ms to years
			if (diff < 1) {
				return d3.rgb("#ff7f0e");
			}
			else {
				return d3.rgb("#1f77b4");
			}
		})	
		.style("stroke", function(d) { 
			var today = new Date();
			var date = new Date(d.date);
			var diff = (today - date) / (1000*60*60*24*365.25) // convert from ms to years
			if (diff < 1) {
				return d3.rgb("#ff7f0e");
			}
			else {
				return d3.rgb("#1f77b4");
			}
		})					
		.on('mouseover', tooltip.show)
      	.on('mouseout', tooltip.hide);
      	
/*	var internalRects = svg.selectAll(".internal")
		.data(internals)
		.enter()
		.append("rect")		
		.attr("class", "internal")
		.attr("width", function(d) { 
			if (typeof d.frequency != "undefined") {
				return 10*Math.sqrt(d.frequency) + 1;
			}
			else {
				return 2;
			}
		})
		.attr("height", function(d) { 
			if (typeof d.frequency != "undefined") {
				return 10*Math.sqrt(d.frequency) + 1;
			}
			else {
				return 2;
			}			
		})		
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
		})			
		.style("fill", function(d) { 
			return d3.rgb("#CCC");
		})		
		.on('mouseover', tooltip.show)
      	.on('mouseout', tooltip.hide)
      	.on('click', function(d) { 
      		var dMin = minimumAttribute(d, "xvalue", d.xvalue),
      			dMax = maximumAttribute(d, "xvalue", d.xvalue),
      			lMin = minimumAttribute(d, "yvalue", d.yvalue),
      			lMax = maximumAttribute(d, "yvalue", d.yvalue);
      		rescale(dMin, dMax, lMin, lMax, xScale, yScale, nodes, links, tips, internals)
      	}); */  	
      	
	d3.select("#reset")
        .on("click", function(d) {
			var dMin = d3.min(xValues),
      			dMax = d3.max(xValues),
      			lMin = d3.min(yValues),
      			lMax = d3.max(yValues);        	
            rescale(dMin, dMax, lMin, lMax, xScale, yScale, nodes, links, tips, internals);
		})  
		
});
