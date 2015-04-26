var ymd_format = d3.time.format("%Y-%m-%d");


var dateValues, earliestDate, dateScale, niceDateScale, counterData;

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

	calcNodeAges(LBItime_window);
	treeplot.selectAll(".link")
		.style("stroke", function(d){return "#ccc";})

	treeplot.selectAll(".tip")
		.style("visibility", tipVisibility)
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
	calcNodeAges(LBItime_window);
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
		.attr("points", branchPoints)
		.style("stroke-width", branchStrokeWidth)
		.style("stroke", branchStrokeColor);				

		d3.selectAll(".tip")
		.transition().duration(500)
		.style("visibility", tipVisibility)
		.style("fill", tipFillColor)
		.style("stroke", tipStrokeColor);
	}
}


function date_init(){

	nodes.forEach(function (d) {d.dateval = new Date(d.date)});
	var dateValues = nodes.filter(function(d) {
		return typeof d.date === 'string';
		}).map(function(d) {
		return new Date(d.date);
	});
	earliestDate = new Date(d3.min(dateValues));
	earliestDate.setDate(earliestDate.getDate() + 180);
	
	var numDateValues = tips.map(function(d) {return d.num_date;})
	var minDate = d3.min(numDateValues);
	var maxDate = d3.max(numDateValues);
	dateDomain = genericDomain.map(function (d){return Math.round(100*(minDate + d*(maxDate - minDate)))/100;});
	dateColorScale.domain(dateDomain);
	dateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true]);	
	
	niceDateScale = d3.time.scale()
		.domain([earliestDate, globalDate])
		.range([5, 205])
		.clamp([true])
		.nice(d3.time.month);
	
	counterData = {}
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

}