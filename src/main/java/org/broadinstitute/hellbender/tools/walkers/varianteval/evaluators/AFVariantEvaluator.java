package org.broadinstitute.hellbender.tools.walkers.varianteval.evaluators;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.hellbender.tools.walkers.varianteval.util.DataPoint;

@Analysis(description = "Calculates the average variant allele fraction (p) and counts the number of each type of variant (Hets/Hom Vars/Hom Refs) over variant sites")

public class AFVariantEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average variant allele fraction over all variant sites", format = "%.8f")
    public double avgVarAlleles= 0.0;
    @DataPoint(description = "Number of called sites over all variant sites;", format = "%d")
    public int totalCalledSites;
    @DataPoint(description = "Number of called heterozygous sites;", format = "%d")
    public int totalHetSites;
    @DataPoint(description = "Number of called homozygous variant sites;", format = "%d")
    public int totalHomVarSites;
    @DataPoint(description = "Number of called homozygous reference sites;", format = "%d")
    public int totalHomRefSites;

    final private double ploidy = 2.0; // assume ploidy of 2
    private double sumVariantAlleles; // this is the allele fraction we're summing over all sites, to be used to calculate the avgVarAlleles


    public AFVariantEvaluator() {
        sumVariantAlleles = 0.0;
        totalCalledSites = 0;
    }

    public int getComparisonOrder() {
        return 1;
    }

    public void update1(VariantContext vc, final ReferenceContext referenceContext, final ReadsContext readsContext, final FeatureContext featureContext) {
        vc.getStart();

        if (vc == null || !vc.isSNP() || (getWalker().ignoreAC0Sites() && vc.isMonomorphicInSamples())) {
            return;
        }

        for (final Genotype genotype : vc.getGenotypes()) {
             // eval array
            if (!genotype.isNoCall()) {
                if (genotype.getPloidy() != ploidy) {

                    throw new UserException.BadInput("This tool only works with ploidy 2");
                }
                // add AF at this site
                this.totalCalledSites += 1;
                int numReferenceAlleles= genotype.countAllele(vc.getReference());
                double varAFHere = (ploidy - numReferenceAlleles)/ploidy;
                this.sumVariantAlleles += varAFHere;

                totalHetSites += numReferenceAlleles == 1 ? 1 : 0;
                totalHomVarSites += numReferenceAlleles == 0 ? 1 : 0;
                totalHomRefSites += numReferenceAlleles == 2 ? 1 : 0;

            }
        }

        if (!vc.hasGenotypes()) {
            // comp  ( sites only thousand genomes )
            this.totalCalledSites += 1;
            this.sumVariantAlleles += vc.getAttributeAsDouble("AF", 0.0);
        }
    }

    @Override
    public void finalizeEvaluation() {
        this.avgVarAlleles = this.totalCalledSites == 0 ? 0 : this.sumVariantAlleles / this.totalCalledSites;
    }
}
