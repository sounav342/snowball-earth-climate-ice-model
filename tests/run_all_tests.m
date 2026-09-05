function results = run_all_tests()
%RUN_ALL_TESTS Run the fast publication-release tests.

root = setup_snowball();
results = runtests(fullfile(root,'tests'),'IncludeSubfolders',true);
assertSuccess(results);
end
