/*
A KBase module: kb_fungalmodeling
This module  build fungal models based on fungal genomes.
*/

module kb_fungalmodeling {


    typedef structure {
        string workspace;
        string genome_ref;
        string template_model;
        int gapfill_model;
        string media_ref;
        string proteintr_ref;
        string translation_policy;
        string custom_model;
        string output_model;
    } fungalmodelbuiltInput;

    typedef structure {
        string report_name;
        string report_ref;
    }fungalmodelbuiltOutput;


    funcdef build_fungal_model(fungalmodelbuiltInput params)
        returns (fungalmodelbuiltOutput output) authentication required;



    typedef structure {
        string workspace;
        string reference_genome;
        string reference_model;
        string genome_ws;
        string model_ws;
    }fungalReferenceModelBuildInput;

    typedef structure {
        string master_template_model_ref;
        string master_template_genome_ref;
    }fungalReferenceModelBuildOutput;

    funcdef build_fungal_template(fungalReferenceModelBuildInput params)
        returns (fungalReferenceModelBuildOutput output) authentication required;

    funcdef build_model_stats(fungalReferenceModelBuildInput params)
        returns (fungalReferenceModelBuildOutput output) authentication required;

    funcdef update_model (fungalReferenceModelBuildInput params)
        returns (fungalReferenceModelBuildOutput output) authentication required;
};
