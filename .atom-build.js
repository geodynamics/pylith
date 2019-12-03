'use babel';

import os  from 'os';

var getBuildDir = function() {
    // :KLUDGE: In order to change to a different build directory, you must restart Atom.
    // It might be better to use the local configuration settings and/or init.coffee.
    if (process.env.PYLITH_BUILDDIR == undefined) {
	atom.notifications.addError('Build directory not set', {
	    detail: "PYLITH_BUILDDIR environment variable not set. Please set PYLITH_BUILDDIR to your pylith build directory and then restart Atom."
	});
    } // if
    
    return process.env.PYLITH_BUILDDIR;
} // getBuildDir


var getMatches = function(output) {
    tags = [
	{
	    type: "Error",
	    regexp: /(.+):(\d+):(\d+): error: (.+)/
	},{
	    type: "Warning",
	    regexp: /(.+):(\d+):(\d+): warning: (.+)/
	}
    ] // tags
	
    var matches = [];
    tags.forEach(tag => {
	output.split('\n').forEach(line => {
	    const match = tag.regexp.exec(line);
	    if (match) {
		matches.push({
		    type: tag.type,
		    file: match[1],
		    line: match[2],
		    col: match[3],
		    message: match[4]
		}); // push
	    } // if
	}); // forEach
    }); // forEach
    return matches;
} // getMatches



var getMakeOptions = function(options) {
    return ["-j"+os.cpus().length].concat(options);
}
	
var createBuild = function(target, options) {
    return {
	name: target,
	cwd: getBuildDir(),
	cmd: "make",
	args: getMakeOptions(options),
	functionMatch: getMatches
    };
} // createBuild

var createDefaultBuild = function(name, options) {
    var all = createBuild("build-all", "");
    all["targets"] = {}
    return all;
} // createDefaultBuild


var addMakeTarget = function(all, name, options) {
    all.targets[name] = createBuild(name, options);
} // addMakeTarget

var builds = function() {
    all = createDefaultBuild("build-all", "all");

    addMakeTarget(all, "build-lib", ["-C libsrc"])
    addMakeTarget(all, "build-module", ["-C modulesrc"])
    addMakeTarget(all, "build-python", ["-C pylith"])
    addMakeTarget(all, "run-libtests", ["check", "-C tests/libtests"])
    addMakeTarget(all, "run-pytests", ["check", "-C tests/pytests"])
    addMakeTarget(all, "run-fullscale", ["check", "-C tests/fullscale"])
    addMakeTarget(all, "run-all-tests", ["check"])
    addMakeTarget(all, "install-all", ["install"])
    addMakeTarget(all, "install-lib", ["install", "-C libsrc"])
    addMakeTarget(all, "install-module", ["install", "-C modulesrc"])
    addMakeTarget(all, "install-python", ["install", "-C pylith"])
    
    return all;
}

module.exports = builds();
