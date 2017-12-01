
package us.kbase.kbfungalmodeling;

import java.util.HashMap;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: fungalmodelbuiltInput</p>
 * 
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "workspace",
    "genome_ref",
    "template_model",
    "gapfill_model",
    "translation_policy",
    "output_model"
})
public class FungalmodelbuiltInput {

    @JsonProperty("workspace")
    private String workspace;
    @JsonProperty("genome_ref")
    private String genomeRef;
    @JsonProperty("template_model")
    private String templateModel;
    @JsonProperty("gapfill_model")
    private Long gapfillModel;
    @JsonProperty("translation_policy")
    private String translationPolicy;
    @JsonProperty("output_model")
    private String outputModel;
    private Map<String, Object> additionalProperties = new HashMap<String, Object>();

    @JsonProperty("workspace")
    public String getWorkspace() {
        return workspace;
    }

    @JsonProperty("workspace")
    public void setWorkspace(String workspace) {
        this.workspace = workspace;
    }

    public FungalmodelbuiltInput withWorkspace(String workspace) {
        this.workspace = workspace;
        return this;
    }

    @JsonProperty("genome_ref")
    public String getGenomeRef() {
        return genomeRef;
    }

    @JsonProperty("genome_ref")
    public void setGenomeRef(String genomeRef) {
        this.genomeRef = genomeRef;
    }

    public FungalmodelbuiltInput withGenomeRef(String genomeRef) {
        this.genomeRef = genomeRef;
        return this;
    }

    @JsonProperty("template_model")
    public String getTemplateModel() {
        return templateModel;
    }

    @JsonProperty("template_model")
    public void setTemplateModel(String templateModel) {
        this.templateModel = templateModel;
    }

    public FungalmodelbuiltInput withTemplateModel(String templateModel) {
        this.templateModel = templateModel;
        return this;
    }

    @JsonProperty("gapfill_model")
    public Long getGapfillModel() {
        return gapfillModel;
    }

    @JsonProperty("gapfill_model")
    public void setGapfillModel(Long gapfillModel) {
        this.gapfillModel = gapfillModel;
    }

    public FungalmodelbuiltInput withGapfillModel(Long gapfillModel) {
        this.gapfillModel = gapfillModel;
        return this;
    }

    @JsonProperty("translation_policy")
    public String getTranslationPolicy() {
        return translationPolicy;
    }

    @JsonProperty("translation_policy")
    public void setTranslationPolicy(String translationPolicy) {
        this.translationPolicy = translationPolicy;
    }

    public FungalmodelbuiltInput withTranslationPolicy(String translationPolicy) {
        this.translationPolicy = translationPolicy;
        return this;
    }

    @JsonProperty("output_model")
    public String getOutputModel() {
        return outputModel;
    }

    @JsonProperty("output_model")
    public void setOutputModel(String outputModel) {
        this.outputModel = outputModel;
    }

    public FungalmodelbuiltInput withOutputModel(String outputModel) {
        this.outputModel = outputModel;
        return this;
    }

    @JsonAnyGetter
    public Map<String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public String toString() {
        return ((((((((((((((("FungalmodelbuiltInput"+" [workspace=")+ workspace)+", genomeRef=")+ genomeRef)+", templateModel=")+ templateModel)+", gapfillModel=")+ gapfillModel)+", translationPolicy=")+ translationPolicy)+", outputModel=")+ outputModel)+", additionalProperties=")+ additionalProperties)+"]");
    }

}
