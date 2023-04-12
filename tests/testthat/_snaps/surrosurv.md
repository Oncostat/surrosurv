# Data

    Code
      nrow(gastadv)
    Output
      [1] 4069
    Code
      names(gastadv)
    Output
      [1] "timeT"    "statusT"  "statusS"  "timeS"    "trialref" "trt"      "id"      

# surrosurv - lite

    {
      "type": "logical",
      "attributes": {
        "dim": {
          "type": "integer",
          "attributes": {},
          "value": [10, 3]
        },
        "dimnames": {
          "type": "list",
          "attributes": {},
          "value": [
            {
              "type": "character",
              "attributes": {},
              "value": ["Clayton unadj", "Clayton adj", "Plackett unadj", "Plackett adj", "Hougaard unadj", "Hougaard adj", "PoissonT", "PoissonI", "PoissonTI", "PoissonTIa"]
            },
            {
              "type": "character",
              "attributes": {},
              "value": ["maxSgrad", "minHev", "minREev"]
            }
          ]
        }
      },
      "value": [false, false, false, false, false, false, false, true, true, true, false, false, false, false, true, true, false, false, true, false, null, false, null, false, null, false, false, null, false, true]
    }

---

    {
      "type": "double",
      "attributes": {
        "dim": {
          "type": "integer",
          "attributes": {},
          "value": [10, 3]
        },
        "dimnames": {
          "type": "list",
          "attributes": {},
          "value": [
            {
              "type": "character",
              "attributes": {},
              "value": ["Clayton unadj", "Clayton adj", "Plackett unadj", "Plackett adj", "Hougaard unadj", "Hougaard adj", "PoissonT", "PoissonI", "PoissonTI", "PoissonTIa"]
            },
            {
              "type": "character",
              "attributes": {},
              "value": ["maxSgrad", "minHev", "minREev"]
            }
          ]
        }
      },
      "value": [10.11028471, 10.11028471, "Inf", "Inf", 31.08446538, 31.08446538, 0.31664406, 0.00185376, 0.0005013, 0.0031577, -0.51101943, -0.51101943, -1, -1, 0.04579231, 0.04579231, -0.00001526, -0.00004578, 0.00019836, -2749747.08578363, "NA", 2.03664888e-09, "NA", -1, "NA", 5.20012419e-13, 2.31282691e-14, "NA", 0, 0.00057605]
    }

---

    {
      "type": "list",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["PoissonTI"]
        },
        "trialSizes": {
          "type": "integer",
          "attributes": {
            "dim": {
              "type": "integer",
              "attributes": {},
              "value": [20]
            },
            "dimnames": {
              "type": "list",
              "attributes": {
                "names": {
                  "type": "character",
                  "attributes": {},
                  "value": [""]
                }
              },
              "value": [
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
                }
              ]
            },
            "class": {
              "type": "character",
              "attributes": {},
              "value": ["table"]
            }
          },
          "value": [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
        },
        "surro.stats": {
          "type": "list",
          "attributes": {
            "dim": {
              "type": "integer",
              "attributes": {},
              "value": [1, 2]
            },
            "dimnames": {
              "type": "list",
              "attributes": {},
              "value": [
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["PoissonTI"]
                },
                {
                  "type": "character",
                  "attributes": {},
                  "value": ["kTau", "R2"]
                }
              ]
            }
          },
          "value": [
            {
              "type": "double",
              "attributes": {},
              "value": [0.53324696]
            },
            {
              "type": "double",
              "attributes": {},
              "value": ["NaN"]
            }
          ]
        }
      },
      "value": [
        {
          "type": "list",
          "attributes": {
            "names": {
              "type": "character",
              "attributes": {},
              "value": ["trtS", "trtT"]
            },
            "class": {
              "type": "character",
              "attributes": {},
              "value": ["data.frame"]
            },
            "row.names": {
              "type": "character",
              "attributes": {},
              "value": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"]
            }
          },
          "value": [
            {
              "type": "double",
              "attributes": {},
              "value": [-0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555, -0.98627555]
            },
            {
              "type": "double",
              "attributes": {},
              "value": [-0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596, -0.54851596]
            }
          ]
        }
      ]
    }

---

    {
      "type": "double",
      "attributes": {
        "names": {
          "type": "character",
          "attributes": {},
          "value": ["Clayton.unadj", "Clayton.adj", "Plackett.unadj", "Plackett.adj", "Hougaard.unadj", "Hougaard.adj", "PoissonT", "PoissonTI", "PoissonTIa"]
        }
      },
      "value": [-1.77702871, -0.97288181, "NA", "NA", -1.77702871, -0.66872289, -0.53821636, "NA", -0.42578614]
    }

