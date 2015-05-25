var legendRectSize = 15;
var legendSpacing = 4;
function makeLegend(){
	
	d3.select("#legend-title").text(function(d){
		if (colorBy == "ep") {
			return "Epitope mutations";
		}
		if (colorBy == "ne") {
			return "Non-epitope mutations";
		}
		if (colorBy == "rb") {
			return "Receptor binding mutations";
		}
		if (colorBy == "lbi") {
			return "Local branching index";
		}
		if (colorBy == "region") {
			return "Region";
		}
		if (colorBy == "genotype") {
			return "Genotype";
		}
		if (colorBy == "date") {
			return "Date";
		}		
		if (colorBy == "dfreq") {
			var tmp_nmonth = Math.round(12*dfreq_dn*time_step);
			var tmp_text = "Freq. change ("+tmp_nmonth+" month";
			if (tmp_nmonth>1){
				tmp_text+='s';
			}
			return tmp_text+')';
		}
});

  legend.attr("height", function(d){return 20*Math.ceil(0.5*colorScale.domain().length);});

  var tmp_leg = legend.selectAll(".legend")
  .data(colorScale.domain())
  .enter().append('g')
  .attr('class', 'legend')
  .attr('transform', function(d, i) {
   var stack = colorScale.domain().length/2;
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
   return d3.rgb(col).toString();
 });

  tmp_leg.append('text')
  .attr('x', legendRectSize + legendSpacing + 5)
  .attr('y', legendRectSize - legendSpacing)
  .text(function(d) {
  	if (typeof d != 'undefined'){
	   return d.toString().replace(/([a-z])([A-Z])/g, '$1 $2').replace(/,/g, ', ');
	}else{ return '';}
 });		
  return tmp_leg;
}

function removeLegend(){
	legend.selectAll('.legend')
  .remove();
}

