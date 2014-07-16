(function ($, ko) {
    'use strict';

    var $container,
        initialized,
        imagesLoaded,
        isotopeOptions,
        indexKey = '__knockout-isotope-index';

    // Modified Isotope methods for gutters in masonry:
    $.Isotope.prototype._getMasonryGutterColumns = function () {
        var gutter = this.options.masonry && this.options.masonry.gutterWidth || 0,
            containerWidth = this.element.width();

        this.masonry.columnWidth = this.options.masonry && this.options.masonry.columnWidth ||
            // or use the size of the first item
            this.$filteredAtoms.outerWidth(true) ||
            // if there's no items, use size of container
            containerWidth;

        this.masonry.columnWidth += gutter;
        this.masonry.cols = Math.floor((containerWidth + gutter) / this.masonry.columnWidth);
        this.masonry.cols = Math.max(this.masonry.cols, 1);
    };

    $.Isotope.prototype._masonryReset = function () {
        // layout-specific props
        this.masonry = {};
        // FIXME shouldn't have to call this again
        this._getMasonryGutterColumns();
        var i = this.masonry.cols;
        this.masonry.colYs = [];
        while (i--) {
            this.masonry.colYs.push(0);
        }
    };

    $.Isotope.prototype._masonryResizeChanged = function () {
        var prevSegments = this.masonry.cols;
        // update cols/rows
        this._getMasonryGutterColumns();
        // return if updated cols/rows is not equal to previous
        return (this.masonry.cols !== prevSegments);
    };

    function afterRender(nodes, data) {
        var $elems = $(nodes).filter(function () { return this.nodeType === 1 });
        $elems.addClass('no-transition');
        $elems.css($container.data('isotope').options.hiddenStyle);
    }

    function afterAdd(node, index, item) {
        if (node.nodeType === 1) { // If this is an html element...
            var $elem = $(node),
                addItem = function () {
                    $elem.data(indexKey, index);
                    $container.isotope('addItems', $elem, function () {
                        $container.isotope().isotope('reLayout');
                        $elem.removeClass('no-transition');
                        $elem.css($container.data('isotope').options.visibleStyle);
                    });
                };

            if (!imagesLoaded)
                $container.imagesLoaded(function () {
                    addItem();
                    imagesLoaded = true;
                });
            else
                addItem();
        }
    }

    function afterMove(node, index, item) {
        if (node.nodeType === 1) // If this is an html element...
            $container.isotope('updateSortData', $(node));
    }

    function beforeRemove(node, index, item) {
        if (node.nodeType === 1) { // If this is an html element...
            var $elem = $(node),
                isotopeInstance = $container.data('isotope');

            isotopeInstance.$allAtoms = isotopeInstance.$allAtoms.not($elem);
            isotopeInstance.$filteredAtoms = isotopeInstance.$filteredAtoms.not($elem);
            $elem.addClass('isotope-hidden').css($container.data('isotope').options.hiddenStyle);
            //isotopeInstance.styleQueue.push({ $el: $elem, style: isotopeInstance.options.hiddenStyle });

            //isotopeInstance._sort(); // dies for some reason... "$.data(...) is undefined" - jquery.isotope.js (line 587)
            $container.isotope('reLayout', function () {
                $elem.remove();
            });
        }
    }

    ko.bindingHandlers.isotope = {
        getSortData: function ($elem) {
            return $elem.data(indexKey);
        },
        makeForeachValueAccessor: function (valueAccessor) {
            return function () {
                var modelValue = valueAccessor(),
                    unwrappedValue = ko.utils.peekObservable(modelValue), // Unwrap without setting a dependency here
                    options = {
                        afterRender: afterRender,
                        afterAdd: afterAdd,
                        afterMove: afterMove,
                        beforeRemove: beforeRemove
                    };

                // If unwrappedValue is the array, pass in the wrapped value on its own
                // The value will be unwrapped and tracked within the template binding
                // (See https://github.com/SteveSanderson/knockout/issues/523)
                if (!unwrappedValue || typeof unwrappedValue.length === "number") {
                    ko.utils.extend(options, {
                        'data': modelValue,
                    });
                    return options;
                }

                // If unwrappedValue.data is the array, preserve all relevant
                // options and unwrap value so we get updates
                ko.utils.unwrapObservable(modelValue);
                ko.utils.extend(options, {
                    'data': unwrappedValue.data,
                });

                return options;
            };
        },
        init: function (element, valueAccessor, allBindingsAccessor, viewModel, bindingContext) {
            var parameters;

            initialized = false;
            isotopeOptions = {};
            $container = $(element);

            parameters = ko.utils.unwrapObservable(valueAccessor());
            if (parameters && typeof parameters == 'object' && !('length' in parameters)) {
                if (parameters.options) {
                    var clientOptions = parameters.options;
                    if (typeof clientOptions !== 'object')
                        throw new Error('options must be an object');
                    $.extend(true, isotopeOptions, clientOptions);
                }
            }

            ko.bindingHandlers.foreach.init(element, ko.bindingHandlers.isotope.makeForeachValueAccessor(valueAccessor), allBindingsAccessor, viewModel, bindingContext);

            return { controlsDescendantBindings: true };
        },
        update: function (element, valueAccessor, allBindingsAccessor, viewModel, bindingContext) {

            ko.bindingHandlers.foreach.update(element, ko.bindingHandlers.isotope.makeForeachValueAccessor(valueAccessor), allBindingsAccessor, viewModel, bindingContext);

            var data = ko.bindingHandlers.isotope.makeForeachValueAccessor(valueAccessor)().foreach;
            ko.utils.unwrapObservable(data);

            if (!initialized) {
                $.extend(true, isotopeOptions, {
                    getSortData: {
                        index: ko.bindingHandlers.isotope.getSortData
                    }
                });

                $container.isotope(isotopeOptions);
            }
            else
                $container.isotope();

            initialized = true;
            return { controlsDescendantBindings: true };
        }
    };
})(window.jQuery, window.ko);