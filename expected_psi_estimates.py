from psi_bayes_estimator import transformed_psi_bayes, curr_psi

# For true psi = .5
print "Current Psi estimate: ", curr_psi(2996, 699, 100, 36, 4)
# If I change overhang to 3, transformed_psi_bayes estimate get much better!
print "Transformed Psi Bayes: ", transformed_psi_bayes(2996, 699, 6304, 300,
                                                       200, 36, 4)
