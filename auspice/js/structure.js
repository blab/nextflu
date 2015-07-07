var myapplett;

function make_structure(){
	console.log('drawing structure');
	var jsmolscript =  "load /data/"+structure+"; cpk off; wireframe off; trace;zoom on;"
					   +"zoom 115;set showhydrogens off; color background white;" 
					   +" select ligand; trace off; spin off; set frank off; "
					   +"set echo bottom left; color echo gray; font echo 14 arial;"
					   +structure_HI_mutations;
//					   +"spacefill 200;color orange;";
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
}

