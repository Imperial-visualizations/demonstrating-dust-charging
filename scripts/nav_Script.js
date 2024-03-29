// Code for navigation functionality in visualisations main screen using Vue

let app = new Vue({

    el: "#app",

    data: {
        // Data required including variables associated with visible sections, script paths and booleans reflecting state of vis
        scrollPos: 0,
        currentTitle: 0,
        currentSection: 0,
        sectionTops: [],
        sectionBottoms: [],
        sectionTitleLong: ["ABR", "OML","MOML", "SOML", "SMOML", "Animation"],
        sectionTitleShort: ["1","2","3","4","5","6"],
        sectionTitle: [],
        hoverPos: '',
        hoverTitle: false,
        mouseX: '',
        n: "",
        journeyHeightOld: "",
        journeyHeightNew: "",
        firstRunDone: false,
        subSection: [false, 1, 1, 1],
    },

    methods: {

        // Function called on scrolling of of left panel to indicate distance scrolled down journey content div
        scrollFunc: function () {
            // function only works once sectionPos has run at least once (in mounted)
            if (app.firstRunDone === true) {
                app.scrollPos = document.querySelectorAll(".journey")[0].scrollTop;
                app.changeTitle();
                app.changeSec();
            }
        },

        handleElement: function (section) {
            // update currentSection variable if user scrolls past the top edge of its corresponding section on left side
            if (app.scrollPos >= app.sectionTops[section - 1] && app.scrollPos < app.sectionBottoms[section - 1]) {
                app.currentTitle = section;
            }
        },

        changeTitle: function () {
            for (let i = 1; i <= app.n; i++) {
                app.handleElement(i)
            }
        },

        changeSec: debounce(function () {
            app.currentSection = app.currentTitle;
        }, 200),

        swapTitles: function (newValue, oldValue) {
            for (let i = 1; i <= app.n; i++) {
                if (i !== newValue) {
                    app.sectionTitle[i - 1] = app.sectionTitleShort[i - 1];
                } else {
                    setTimeout(function () {
                        app.sectionTitle[i - 1] = app.sectionTitleLong[i - 1];
                    }, 20);
                    setTimeout(function () {
                        app.$forceUpdate();
                    }, 100);
                }
            }
        },

        // Function called every x seconds to check if section div sizes have changed and recalculate scroll positions if so
        // Div sizes may change if window re-sized or if a subsection is expanded/collapsed
        sectionPos: function () {
            this.$nextTick(function () {
                let overallTop = document.querySelectorAll("#sc1")[0].offsetTop;
                for (let i = 1; i <= app.n; i++) {
                    if (i < app.n) {
                        app.sectionTops[i - 1] = (document.querySelectorAll("#" + "sc" + i)[0].offsetTop - overallTop);
                        app.sectionBottoms[i - 1] = (app.sectionTops[i - 1] + document.querySelectorAll("#" + "sc" + i)[0].offsetHeight);
                    } else {
                        app.sectionTops[i - 1] = (document.querySelectorAll("#" + "sc" + i)[0].offsetTop - overallTop);
                        app.sectionBottoms[i - 1] = (app.sectionTops[i - 1] + document.querySelectorAll("#" + "sc" + i)[0].offsetHeight - document.querySelectorAll(".journey")[0].offsetHeight);
                    }
                }
                app.firstRunDone = true;
                app.scrollFunc();
            })
        },

        // Function activated when button in nav/progress bar clicked to scroll automatically to relevant section
        scrollTo: function (event) {
            document.querySelectorAll("#" + "ph" + event.currentTarget.dataset.no)[0].scrollIntoView({behavior: "smooth"});
        },

        // Same as above but for subsections
        subScrollTo: function (section) {
            if (app.currentSection !== section) {
                let scrollTarget = document.querySelectorAll("#ph" + section)[0];
                scrollTarget.scrollIntoView({behavior: "smooth"});
            }
        },

        // Updates number of title being hovered over in nav/progress bar in data
        hoverPosUpdate: function (event) {
            app.hoverPos = parseFloat(event.currentTarget.dataset.no)
        },

        // Updates if and what title show when hovering over nav/progress bar
        selectHover: function () {
            if (app.currentTitle !== app.hoverPos) {
                app.hoverTitle = app.sectionTitleLong[app.hoverPos - 1]
            } else {
                app.hoverTitle = false
            }
        },

        // Updates x-position of mouse in data
        updateMouseX: function (event) {
            // pass event object, bound to mouse move with update
            app.mouseX = event.clientX - 15;
        },

    },

    watch: {

        // Updates current section title to display in full in nav/progress bar whilst minimising other section titles
        currentTitle: function (newValue, oldValue) {
            app.swapTitles(newValue, oldValue)
        },
    },

    mounted() {
        // $nextTick ensures initial functions only run once Vue is initialised sufficiently
        this.$nextTick(function () {
                // makes n equal to total number of sections
                app.n = document.querySelectorAll(".section-container").length;
                // calculates initial div section positions in journey with respect to the top
                app.sectionPos();
                // checks if journey div height changes every x seconds
                // if it does change, re-runs sectionPos to calculate section div positions
                app.journeyHeightOld = document.querySelectorAll(".journey")[0].scrollHeight;
                window.setInterval(() => {
                    app.journeyHeightNew = document.querySelectorAll(".journey")[0].scrollHeight;
                    if (app.journeyHeightOld !== app.journeyHeightNew) {
                        app.journeyHeightOld = app.journeyHeightNew;
                        this.sectionPos();
                    }
                }, 2000);
            }
        )
    },
});