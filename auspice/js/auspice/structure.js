//var options = {
//  width: 500,
//  height: 500,
//  antialias: true,
//  quality : 'medium'
//};
//
//function make_structure(){
//	// insert the viewer under the Dom element with id 'gl'.
//	var parent = document.getElementById('HA_struct')
//	var viewer = pv.Viewer(parent, options);
//
//	function setColorForAtom(go, atom, color) {
//	    var view = go.structure().createEmptyView();
//	    view.addAtom(atom);
//	    go.colorBy(pv.color.uniform(color), view);
//	}
//
//	// variable to store the previously picked atom. Required for resetting the color
//	// whenever the mouse moves.
//	var prevPicked = null;
//	// add mouse move event listener to the div element containing the viewer. Whenever
//	// the mouse moves, use viewer.pick() to get the current atom under the cursor.
//	parent.addEventListener('mousemove', function(event) {
//	    var rect = viewer.boundingClientRect();
//	    var picked = viewer.pick({ x : event.clientX - rect.left,
//	                               y : event.clientY - rect.top });
//	    if (prevPicked !== null && picked !== null &&
//	        picked.target() === prevPicked.atom) {
//	      return;
//	    }
//	    if (prevPicked !== null) {
//	      // reset color of previously picked atom.
//	      setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
//	    }
//	    if (picked !== null) {
//	      var atom = picked.target();
//	      document.getElementById('residue_name').innerHTML = atom.qualifiedName();
//	      // get RGBA color and store in the color array, so we know what it was
//	      // before changing it to the highlight color.
//	      var color = [0,0,0,0];
//	      picked.node().getColorForAtom(atom, color);
//	      prevPicked = { atom : atom, color : color, node : picked.node() };
//
//	      setColorForAtom(picked.node(), atom, 'red');
//	    } else {
//	      document.getElementById('residue_name').innerHTML = '&nbsp;';
//	      prevPicked = null;
//	    }
//	    viewer.requestRedraw();
//	});
//	pv.io.fetchPdb("/data/"+structure, function(structure) {
//	  // display the protein as cartoon, coloring the secondary structure
//	  // elements in a rainbow gradient.
//	  viewer.cartoon('protein', structure); //, { color : color.ssSuccession() });
//	  viewer.centerOn(structure);
//	  viewer.on('viewerReady', function() {
//    	  var go = viewer.cartoon('structure', structure);
//      	// adjust center of view and zoom such that all structures can be seen.
//      	viewer.autoZoom();
//	   });o
//	});
//}
//
//var myapplett;
//
function make_structure(){
	console.log('drawing structure');
	var jsmolscript =  "load /data/"+structure+"; cpk off; wireframe off; cartoon ONLY; trace;zoom on;"
					   +"zoom 115;set showhydrogens off; color background white;"
					   +" select ligand; trace off; spin off; set frank off; "
					   +"set echo bottom left; color echo gray; font echo 14 arial;"
					   +structure_HI_mutations
					   +" select (chain==C); color [xFFFFFF]; select (chain==D); color [xFFFFFF];"
					   +" select (chain==E); color [xFFFFFF]; select (chain==F); color [xFFFFFF];";
	console.log(jsmolscript);
	Info = {
		width: 500,
		height: 500,
		debug: false,
		j2sPath: "/js/j2s",
		color: "white",
		disableJ2SLoadMonitor: true,
		disableInitialConsole: true,
		addSelectionOptions: false,
		use: "HTML5",
		readyFunction: null,
		script:	jsmolscript}

	myapplett = $("#HA_struct").html(Jmol.getAppletHtml("jmolApplet0",Info));
	var structCaption = document.getElementById('struct_caption');
	struct_caption.innerHTML='JSmol rendering of <a target="_blank" href="http://www.rcsb.org/pdb/explore/explore.do?structureId='+structure.substring(0,4)+'">'+structure.substring(0,4)+'</a>';
}

