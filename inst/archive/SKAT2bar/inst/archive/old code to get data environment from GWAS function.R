#noGform= noCallPattern(formula,call.pattern=c(as.name(".eigenG"),as.name(".G")))
#noGformGxE=noCallPattern(GxE.formula)
#trueTerms=c(findTrueTerms(noGform),findTrueTerms(noGformGxE))
#trueTerms.form=formula
#RHSForm(trueTerms.form)=asCall(paste0(trueTerms,collapse="+"))
#combinedFormMto1=replaceTerm(trueTerms.form,movingTerm,1)

##get data environment for all the variables not in the moving window. Therefore, all the terms will be na.omit at the same time except for the moving term, which will be treated separately.
#   if (is.null(data)) {
#     fr=model.frame(combinedFormMto1)
#   } else {
#     fr=model.frame(combinedFormMto1,data=data)
#   }
#   names(fr)=gsub(".*\\$(.*)$","\\1",names(fr))
#   denv=list2env(fr)
#   whichNa=as.integer(attr(fr,"na.action"))
#   n0=nrow(fr)+length(whichNa)
#   rm("fr")
#   formula=removeDollar(formula)
#   environment(formula)=denv
#   rm("combinedFormMto1")
#   rm("trueTerms")
#   rm("trueTerms.form")