/*
A KBase module: kb_fungalmodeling
This sample module contains one small method - filter_contigs.
*/

module kb_fungalmodeling {


    typedef structure {
        string  workspace;
        string  genome_ref;
        string  template_model;
        string  translation_policy;
        string  output_model;
    } fungalmodelbuiltInput;

    typedef structure {
        string report_name;
        string report_ref;


    }fungalmodelbuiltOutput;


    funcdef build_fungal_model(fungalmodelbuiltInput params)
        returns (fungalmodelbuiltOutput output) authentication required;
};
