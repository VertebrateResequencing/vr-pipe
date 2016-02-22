(function(factory) {"use strict";
	if (typeof require === "function" && typeof exports === "object" && typeof module === "object") {
		// CommonJS or Node: hard-coded dependency on "knockout"
		factory(require("knockout"), exports);
	} else if (typeof define === "function" && define["amd"]) {
		// AMD anonymous module with hard-coded dependency on "knockout"
		define(["knockout", "exports"], factory);
	} else {
		// <script> tag: use the global `ko` object
		factory(ko, {});
	}
}(function(ko, exports) {"use strict";
	ko.extenders.progressivefilter = function(target, args) {
		var requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame || window.webkitRequestAnimationFrame || window.msRequestAnimationFrame || function(callback) { setTimeout(callback, 4); },
			currentCount = 0,
			props = {};

		args = args || {};
		target.progressivefilter = props;

		props.unfilteredCollection = [];
		props.unfilteredCollectionIndex = 0;
		props.isFiltering = ko.observable(false);
		props.filterFunction = args.filterFunction;
		props.batchSize = Math.max(parseInt(args.batchSize, 10), 1);

		props.add = args.addFunction || function(item) { target.peek().push(item); };
        props.clear = args.clearFunction || function() { target([]); };

		target.isFiltered = function(item) {
			return !props.filterFunction || props.filterFunction(item);
		};

		target.filter = function(unfilteredCollection) {
			var filteredCollection = [],
				i;
			for (i = 0; i < unfilteredCollection.length; i++) {
				if (target.isFiltered(unfilteredCollection[i])) {
					filteredCollection.push(unfilteredCollection[i]);
				}
			}
			props.clear();
			target(filteredCollection);
		};

		target.filterProgressive = function(unfilteredCollection) {
			props.unfilteredCollection = unfilteredCollection.slice(0);
			props.unfilteredCollectionIndex = 0;
			currentCount = 0;
			props.clear();
			if (!props.isFiltering.peek()) {
				props.isFiltering(true);
				requestAnimationFrame(doFilter);
			}
		};

		function doFilter() {
			var item;

			for (props.unfilteredCollectionIndex; props.unfilteredCollectionIndex < props.unfilteredCollection.length; props.unfilteredCollectionIndex++) {
				item = props.unfilteredCollection[props.unfilteredCollectionIndex];
				if (item && target.isFiltered(item)) {
					props.add(item);
					break;
				}
			}

			currentCount++;
			props.unfilteredCollectionIndex++;

			if (props.unfilteredCollectionIndex < props.unfilteredCollection.length) {
				if (currentCount >= props.batchSize) {
					target.valueHasMutated();
					currentCount = 0;
					requestAnimationFrame(doFilter);
				}
				else {
					currentCount++;
					doFilter();
				}
				return;
			}
			else {
                target.valueHasMutated();
                currentCount = 0;
				props.unfilteredCollectionIndex = 0;
				props.isFiltering(false);
			}
		}
	};
}));