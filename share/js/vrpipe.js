// we have a custom ko binding to handle our loading feedback
// (from http://blog.greatrexpectations.com/2012/06/17/loading-placeholders-using-knockout-js/)
ko.bindingHandlers.loadingWhen = {
    init: function (element) {
        var $element = $(element),
            currentPosition = $element.css("position")
            $loader = $("<div>").addClass("loader").hide();
        
        //add the loader
        $element.append($loader);
        
        //make sure that we can absolutely position the loader against the original element
        if (currentPosition == "auto" || currentPosition == "static")
            $element.css("position", "relative");

        //center the loader
        $loader.css({
            position: "absolute",
            top: "20px",
            left: "50%",
            "margin-left": -($loader.width() / 2) + "px"
        });
    },
    update: function (element, valueAccessor) {
        var isLoading = ko.utils.unwrapObservable(valueAccessor()),
            $element = $(element),
            $childrenToHide = $element.children(":not(div.loader)"),
            $loader = $element.find("div.loader");

        if (isLoading) {
            $childrenToHide.css("visibility", "hidden").attr("disabled", "disabled");
            $loader.fadeIn("slow");
        }
        else {
            $loader.fadeOut("fast");
            $childrenToHide.css("visibility", "visible").removeAttr("disabled");
        }
    }
};

// we can have sortable table columns
// (from http://stackoverflow.com/a/16964843/675083)
ko.bindingHandlers.sort = {
    init: function(element, valueAccessor, allBindingsAccessor, viewModel, bindingContext) {
        var asc = false;
        element.style.cursor = 'pointer';
        
        element.onclick = function(){
            var value = valueAccessor();
            var prop = value.prop;
            var data = value.arr;
            
            asc = !asc;
            
            data.sort(function(left, right){
                var rec1 = left;
                var rec2 = right;
                
                if(!asc) {
                    rec1 = right;
                    rec2 = left;
                }
                
                var props = prop.split('.');
                for(var i in props){
                    var propName = props[i];
                    var parenIndex = propName.indexOf('()');
                    if(parenIndex > 0){
                        propName = propName.substring(0, parenIndex);
                        rec1 = rec1[propName]();
                        rec2 = rec2[propName]();
                    } else {
                        rec1 = rec1[propName];
                        rec2 = rec2[propName];
                    }
                }
                
                return rec1 == rec2 ? 0 : rec1 < rec2 ? -1 : 1;
            });
        };
    }
};

// load knockout templates from an external file
// from http://stackoverflow.com/questions/9480697/using-an-external-template-in-knockoutjs/22902067#22902067
// use like:
// <div id="nav">
//     <script type="text/html" src="/login_template.html" id="login-template"></script>
//     <div data-bind="template: { name: 'login-template' }"></div>
// </div>
// var livm = new ko.LogInViewModel();
// loadExternalKnockoutTemplates(function() {
//     ko.applyBindings(livm, $('#nav')[0]);
// });
function loadExternalKnockoutTemplates(callback) {
    var sel = 'script[src][type="text/html"]:not([loaded])';
    $toload = $(sel);
    function oncomplete() {
        this.attr('loaded', true);
        var $not_loaded = $(sel);
        if(!$not_loaded.length) {
            callback();
        }
    }
    ko.utils.arrayForEach($toload, function(elem) {
        var $elem = $(elem);
        $elem.load($elem.attr('src'), oncomplete.bind($elem));

    });
}

// a knockout view model for logging in, which can be 'inherited' in another
// viewmodel (for eg. your navbar) by doing:
// ko.LogInViewModel.call(self);
// Or to just use it directly:
// var livm = new ko.LogInViewModel(); ko.applyBindings(livm, $('#nav')[0]);
// If the user successfully authenticates, the username() observable will
// contain their username.
// You must have html elements with ids loginUsername, loginPassword,
// loginSubmit and loginDropdown on your page for this to work. You should
// show the loginError() observable somewhere as well. Best to just use the
// login_template.html (see comments for loadExternalKnockoutTemplates).
(function (ko, undefined) {
    ko.LogInViewModel = function () {
        var self = this;
        self.loading = ko.observableArray();
        self.errors = ko.observableArray();
        self.username = ko.observable();
        self.loginError = ko.observable();
        
        self.getUserName = function() {
            // ajax call to refresh the idle session timeout and populate
            // username if user has already logged in
            vrpipeRestMethod('authenticated_user', 'n/a', {}, self.loading, self.errors, function(data) {
                self.username(data['user']);
            });
        }
        
        // to set username immediately if logged in, and to keep the session
        // alive if we are already logged in, first we call getUserName() on
        // page load, then we call it every hour if there has been any activity
        // on the page in the last hour
        self.getUserName();
        self.intervalMilliSeconds = 3600000;
        self.hadActivity = false;
        $(document.body).bind("mousemove keypress", function(e) {
            self.hadActivity = true;
        });
        self.maintainSession = function() {
            if (self.hadActivity) {
                self.getUserName();
            }
            self.hadActivity = false;
        }
        setInterval(self.maintainSession, self.intervalMilliSeconds);
        
        self.login = function() {
            var user = $("#loginUsername").val();
            var password = $("#loginPassword").val();
            if (user && password) {
                $('#loginSubmit').prop('disabled', true);
                vrpipeRestMethod('authenticate', 'n/a', { user: user, password: password }, self.loading, self.errors, function(data) {
                    $('#loginSubmit').prop('disabled', false);
                    if (data['user']) {
                        $('#loginDropdown').dropdown('toggle');
                        self.username(data['user']);
                        window.location.reload();
                    }
                    else {
                        self.loginError('Failed; try again?');
                    }
                });
            }
        }
        
        self.logout = function() {
            vrpipeRestMethod('deauthenticate_user', 'n/a', { }, self.loading, self.errors, function(data) {
                self.loginError(undefined);
                self.username(undefined);
                window.location.reload();
            });
        }
    };
}(ko));

// ajax function to send/get json queries
var ajax = function(uri, method, loading, error, data) {
    var loadingTimeout;
    var request = {
        url: uri,
        type: method,
        contentType: "application/json",
        accepts: "application/json",
        cache: false,
        dataType: 'json',
        data: JSON.stringify(data),
        beforeSend: function (xhr) {
            // say that we're loading, but only if for more than
            // half a second
            if (loadingTimeout) { clearTimeout(loadingTimeout); }
            loadingTimeout = setTimeout(function() {
                loading.push(true);
            }, 500);
        },
        error: function(jqXHR, textStatus, errorThrown) {
            if (textStatus) {
                var msg = "Ajax " + textStatus + " " + jqXHR.status + " connecting to " + uri;
                if (errorThrown) {
                    msg = msg + ": " + errorThrown;
                }
                error.push(msg);
            }
        },
        complete: function(jqXHR, textStatus) {
            if (loadingTimeout) { clearTimeout(loadingTimeout); }
            loading.pop();
        }
    };
    return $.ajax(request);
}

// function to call a VRPipe Rest method and do something with the
// json response, and its helper function to deal with errors
var vrpipeRestDataErrorParser = function(data, error) {
    if (data == null) {
        error.push('no data');
        return 1;
    }
    if (data.hasOwnProperty('errors')) {
        for (var i = 0; i < data.errors.length; i++) {
            error.push(data.errors[i]);
        }
        return 1;
    }
    return 0;
}
var vrpipeRestMethod = function(category, method, args, loading, error, donefunction, other) {
    var promise = ajax('/rest/' + category + '/' + method, 'POST', loading, error, args);
    
    if (donefunction) {
        promise.done(function(data) {
            if (! vrpipeRestDataErrorParser(data, error)) {
                donefunction(data, other);
            }
        });
    }
    else {
        return promise;
    }
}

var shortenString = function(string) {
    var shortened = string;
    if (string.length > 25) {
        shortened = string.substring(0, 5) + '[...]' + string.slice(-15);
    }
    return shortened;
}

var wbr = function(string) {
    if (! string) { return '' }
    var wbrd = '' + string;
    wbrd = wbrd.replace(/([_\-\W])/g, "$1<wbr>");
    wbrd = wbrd.replace(/( <wbr>)/g, " ");
    return wbrd;
}

var zeropad = function(num) {
    return (num < 10) ? ("0" + num) : num;
}

// convert epoch seconds to friendly date
var epochToDate = function(epoch) {
    var date = new Date(epoch * 1000);
    return date.getFullYear() + "/" + zeropad(date.getMonth() + 1) + "/" + zeropad(date.getDate()) + " " + zeropad(date.getHours()) + ":" + zeropad(date.getMinutes()) + ":" + zeropad(date.getSeconds());
}

// find out if 2 arrays are equal, regardless of order
function arraysEqual(a, b) {
    if (a === b) return true;
    if (a == null || b == null) return false;
    if (a.length != b.length) return false;

    a.sort();
    b.sort();

    for (var i = 0; i < a.length; ++i) {
        if (a[i] !== b[i]) return false;
    }
    return true;
}

// split a comma separated string of numbers in to an array of ints
function splitToInts(str) {
    return $.map(str.split(','), function(value){return parseInt(value, 10);});
}

// wait until a given observeable array has length, then run the given function
var runWhenPopulated = function(observable, runFunction) {
    if (observable().length) {
        runFunction();
    }
    else {
        var subscription;
        subscription = observable.subscribe(function(newValue) {
            if (newValue.length) {
                subscription.dispose();
                runFunction();
            }
        });
    }
}

// get a stacktrace
var stackTrace = function() {
    var err = new Error();
    console.log(err.stack);
}